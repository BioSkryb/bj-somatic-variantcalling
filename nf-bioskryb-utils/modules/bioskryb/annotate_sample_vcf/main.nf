nextflow.enable.dsl=2
params.timestamp = ""

// For one sample: extract a per-sample VCF from the group VCF, then annotate it
// with INFO fields derived from:
//   (a) vcf_annotation_table_${group}.tsv — per-variant filter provenance (29 fields)
//   (b) pileup_focal_variants_${group}.tsv — per-sample pileup metrics (24 fields, SMPL_PILEUP_ prefix)
// Both tables are joined by awk into a single combined TSV so that bcftools annotate
// runs only once (53 INFO fields total).
// One Nextflow task per sample for maximum parallelism.
process ANNOTATE_SAMPLE_VCF {
    tag "${group}:${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), val(sample_name),
          path(group_vcf), path(tbi),
          path(annot_table),
          path(pileup_focal)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), val(sample_name),
          path("${sample_name}_somatic_annotated.vcf.gz"),
          path("${sample_name}_somatic_annotated.vcf.gz.tbi"),
          emit: annotated_vcf

    script:
    """
    set -euo pipefail

    # ── A: Join both annotation sources into a single combined TSV ───────────────
    # Built first so combined.tsv.gz can serve as a targets file for bcftools view.
    # File 1 (pileup_focal): cols 1=SampleId 2=VariantId 3=CHROM 4=POS 5=REF 6=ALT 7..=metrics
    #   → load rows for this sample into a lookup keyed by VariantId
    # File 2 (annot_table): col 1=VariantId (CHROM_POS_REF_ALT format) 2..=filter fields
    #   → parse coords from VariantId, emit merged row (NA for missing pileup entries)
    # VariantId parsing: last field=ALT, n-1=REF, n-2=POS, everything before=CHROM
    # (handles non-standard names like chrUn_gl000220).
    awk -F'\\t' -v smp="${sample_name}" '
    FILENAME == ARGV[1] && FNR == 1 {
        n_pileup = NF - 6
        for (i = 7; i <= NF; i++) pileup_hdr[i-6] = "SMPL_PILEUP_" \$i
        next
    }
    FILENAME == ARGV[1] && \$1 == smp && \$6 != "REF" {
        vid = \$2
        val = \$7
        for (i = 8; i <= NF; i++) val = val "\\t" \$i
        pileup_alt[vid] = val
        next
    }
    FILENAME == ARGV[1] && \$1 == smp && \$6 == "REF" {
        vid = \$2
        # Cols 7-10 are allele-specific fragment counts (REF allele reads in a REF row).
        # For 0/0 sites ALT support must be 0; keep col 11+ (quality metrics, depth, filters).
        val = "0\\t0\\t0\\t0"
        for (i = 11; i <= NF; i++) val = val "\\t" \$i
        pileup_ref[vid] = val
        next
    }
    FILENAME == ARGV[2] && FNR == 1 {
        printf "CHROM\\tPOS\\tREF\\tALT"
        for (i = 2; i <= NF; i++) printf "\\t%s", \$i
        for (i = 1; i <= n_pileup; i++) printf "\\t%s", pileup_hdr[i]
        print ""
        # build NA placeholder for variants absent in pileup for this sample
        na_str = "NA"
        for (i = 2; i <= n_pileup; i++) na_str = na_str "\\tNA"
        next
    }
    FILENAME == ARGV[2] {
        n = split(\$1, a, "_")
        alt = a[n]; ref = a[n-1]; pos = a[n-2]
        chrom = a[1]; for (j = 2; j <= n-3; j++) chrom = chrom "_" a[j]
        printf "%s\\t%s\\t%s\\t%s", chrom, pos, ref, alt
        for (i = 2; i <= NF; i++) printf "\\t%s", \$i
        if      (\$1 in pileup_alt) printf "\\t%s", pileup_alt[\$1]
        else if (\$1 in pileup_ref) printf "\\t%s", pileup_ref[\$1]
        else                         printf "\\t%s", na_str
        print ""
    }
    ' ${pileup_focal} ${annot_table} \
    | sort -k1,1 -k2,2n \
    | bgzip -c > combined.tsv.gz
    tabix -s1 -b2 -e2 -S1 combined.tsv.gz

    # ── B: Extract per-sample VCF with exact CHROM_POS_REF_ALT matching ──────────
    # Step 1: regions.bed (CHROM, start-1, end) for bcftools -R broad prefilter.
    zcat combined.tsv.gz | awk 'NR>1 {print \$1"\\t"(\$2-1)"\\t"\$2}' \
        | sort -k1,1 -k2,2n | uniq > regions.bed

    # Step 2: exact CHROM_POS_REF_ALT lookup (cols 1-4 of combined.tsv.gz).
    zcat combined.tsv.gz | awk 'NR>1 {print \$1"_"\$2"_"\$3"_"\$4}' > focal_variants.txt

    # Step 3: broad region extract (all genotypes including 0/0) → normalise ./. → 0/0
    #         → exact-match filter → per-sample VCF with all 668 focal positions.
    bcftools view -s ${sample_name} -R regions.bed -Ov ${group_vcf} \
        | bcftools +setGT -- -t ./. -n 0/0 \
        > variants_tmp.vcf
    grep "^#"  variants_tmp.vcf > header.vcf
    grep -v "^#" variants_tmp.vcf > body.vcf
    awk -F'\\t' -v OFS='\\t' \
        'NR==FNR { a[\$0]; next } { if ((\$1"_"\$2"_"\$4"_"\$5) in a) print }' \
        focal_variants.txt body.vcf > chosen_body.vcf
    cat header.vcf chosen_body.vcf | bcftools view -Oz -o ${sample_name}_raw.vcf.gz
    bcftools index -t ${sample_name}_raw.vcf.gz

    # ── C: Build ##INFO header lines and -c column mapping ───────────────────────
    # Read the header once; pipefail is disabled for this one pipeline because
    # head -1 closes the pipe early, causing zcat to exit with SIGPIPE (141).
    set +o pipefail
    HEADER_FIELDS=\$(zcat combined.tsv.gz | head -1 | cut -f5- | tr '\\t' '\\n')
    set -o pipefail

    echo "\${HEADER_FIELDS}" | awk '
        /^SMPL_PILEUP_/ { print "##INFO=<ID=" \$1 ",Number=1,Type=String,Description=\\"Per-sample pileup metric: " \$1 "\\">"; next }
                        { print "##INFO=<ID=" \$1 ",Number=1,Type=String,Description=\\"Somatic filter provenance: " \$1 "\\">"}
    ' > combined_hdr.txt

    ALL_COLS=\$(echo "\${HEADER_FIELDS}" | awk '{printf ",INFO/%s", \$1}' | sed 's/^,//')

    # ── D: Single bcftools annotate pass (all 53 INFO fields) ────────────────────
    bcftools annotate \
        -a combined.tsv.gz \
        -h combined_hdr.txt \
        -c CHROM,POS,REF,ALT,\${ALL_COLS} \
        -Oz -o ${sample_name}_somatic_annotated.vcf.gz \
        ${sample_name}_raw.vcf.gz
    bcftools index -t ${sample_name}_somatic_annotated.vcf.gz

    echo "Done: \$(bcftools view -H ${sample_name}_somatic_annotated.vcf.gz | wc -l) variants annotated for ${sample_name}"
    """
}
