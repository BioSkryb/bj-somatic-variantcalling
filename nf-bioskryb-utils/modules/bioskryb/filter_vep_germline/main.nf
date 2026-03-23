nextflow.enable.dsl=2
params.timestamp = ""

// Filter VEP-annotated VCF: optionally exclude variants only in dbSNP (Existing_variation);
// always apply AF/MAX_AF <= max_af. When filter_by_existing_variation is true: keep only if no IDs, or ≥1 non-dbSNP ID.
// Input VCF is sorted and indexed internally before filtering.
// Outputs:
//   chosen_variants: CHROM_POS_REF_ALT (one per line) of variants that passed all filters
//   filter_provenance: TSV table of all variants with FILTER_STATUS (PASS/REMOVED) and FILTER_REASON
process FILTER_VEP_GERMLINE {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), path(vep_vcf)
    val(max_af)
    val(filter_by_existing_variation)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), path("vep_sorted_${group}.vcf.gz"), path("vep_sorted_${group}.vcf.gz.tbi"), emit: sorted_vep_vcf
    tuple val(group), path("chosen_variants_postgermlinefilter_${group}.txt"), emit: chosen_variants
    tuple val(group), path("vep_filter_provenance_${group}.tsv"), emit: filter_provenance

    script:
    def vcf = vep_vcf[0]
    def filter_ev = filter_by_existing_variation ? 1 : 0
    // bcftools -i/-e does not see CSQ-extracted AF when using -f (format output), so do AF + optional Existing_variation in awk.
    // -d: decompose multi-allelic so each row is a single ALT. -s worst: one CSQ row per variant (worst consequence).
    // Cols in split_vep_raw.tsv: CHROM, POS, REF, ALT, Existing_variation, AF, MAX_AF
    """
    bcftools sort -Oz -o vep_sorted_${group}.vcf.gz ${vcf}
    bcftools index -t vep_sorted_${group}.vcf.gz

    echo -e "Extracting VEP fields ..."
    bcftools +split-vep vep_sorted_${group}.vcf.gz \\
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%Existing_variation\\t%AF\\t%MAX_AF\\n' \\
        -d -s worst > split_vep_raw.tsv

    echo -e "Filtering VEP VCF: filter_by_existing_variation=${filter_by_existing_variation}; exclude AF, MAX_AF > ${max_af} ..."

    echo -e "CHROM\\tPOS\\tREF\\tALT\\tExisting_variation\\tAF\\tMAX_AF\\tFILTER_STATUS\\tFILTER_REASON" \\
        > vep_filter_provenance_${group}.tsv

    awk -v max_af=${max_af} -v filter_ev=${filter_ev} \\
        -v prov="vep_filter_provenance_${group}.tsv" \\
        -v chosen="chosen_variants_postgermlinefilter_${group}.txt" \\
        'BEGIN{FS="\\t"; OFS="\\t"}
      NF>=7 {
        keep_ev=1; keep_af=1; reason="";

        # --- Existing_variation filter ---
        if (filter_ev) {
          ev=\$5; keep_ev=0;
          if (ev=="" || ev==".") keep_ev=1;
          else {
            gsub(/&/, ",", ev); n=split(ev,a,",");
            for(i=1;i<=n;i++) if (a[i]!="" && a[i]!~/^rs[0-9]+\$/) { keep_ev=1; break }
          }
        }

        # --- AF / MAX_AF filter ---
        if (!((\$6=="" || \$6=="." || \$6+0<=max_af) && (\$7=="" || \$7=="." || \$7+0<=max_af))) keep_af=0;

        # --- Determine status and reason ---
        if (!keep_ev) reason="existing_variation_filter";
        if (!keep_af) reason=(reason=="" ? "AF_filter" : reason";AF_filter");
        status=(keep_ev && keep_af) ? "PASS" : "REMOVED";
        if (status=="PASS") reason=".";

        # --- Write provenance row ---
        print \$1,\$2,\$3,\$4,\$5,\$6,\$7,status,reason >> prov;

        # --- Write to chosen_variants if PASS ---
        if (status=="PASS") print \$1"_"\$2"_"\$3"_"\$4 >> chosen;
      }' split_vep_raw.tsv

    sort -u chosen_variants_postgermlinefilter_${group}.txt \\
        -o chosen_variants_postgermlinefilter_${group}.txt

    echo -e "Wrote \$(wc -l < chosen_variants_postgermlinefilter_${group}.txt) variants to chosen_variants_postgermlinefilter_${group}.txt"
    echo -e "Wrote \$(wc -l < vep_filter_provenance_${group}.tsv) lines (incl. header) to vep_filter_provenance_${group}.tsv"
    """
}
