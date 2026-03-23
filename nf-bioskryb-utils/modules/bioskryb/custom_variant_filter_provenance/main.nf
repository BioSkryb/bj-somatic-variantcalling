nextflow.enable.dsl=2
params.timestamp = ""

// Build a master variant × filter-status table by joining all per-stage provenance outputs.
// Also generates a filter tracking table showing how many variants are retained at each stage.
// Join order: bulk → binom → vep → pileup → sequoia.
// All columns from every source table are kept, prefixed with the source name.
// Pileup is sample-level: numeric metrics are averaged, filter cols become PassCounts,
// SampleId becomes NumSamples, and per-variant-constant fields take their first value.
process CUSTOM_VARIANT_FILTER_PROVENANCE {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/CUSTOM_VARIANT_FILTER_PROVENANCE/", enabled: "$enable_publish"

    input:
    tuple val(group), path(merged_vcf_variants), path(bulk_prov), path(binom_tsv), path(vep_prov), path(pileup_tables), path(sequoia_filt), path(tab_nvnr_files)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), path("variant_master_filter_table_${group}.tsv"),        emit: master_table
    tuple val(group), path("variant_filter_tracking_${group}.tsv"),            emit: filter_tracking
    tuple val(group), path("variant_filter_report_${group}.md"),               emit: filter_report
    tuple val(group), path("variant_filter_plot_${group}.pdf"),                emit: filter_plot
    tuple val(group), path("binom_filter_plots_${group}.pdf"),                 emit: binom_plots
    tuple val(group), path("vaf_hexbin_plots_${group}.pdf"),                   emit: vaf_hexbin_plots
    tuple val(group), path("pileup_metric_plots_${group}.pdf"),                emit: pileup_metric_plots
    tuple val(group), path("pileup_bppos_plots_${group}.pdf"),                 emit: pileup_bppos_plots
    tuple val(group), path("Patient_filter_report_combined_${group}.pdf"),     emit: combined_report
    tuple val(group), path("vcf_annotation_table_${group}.tsv"),               emit: vcf_annotation_table
    tuple val(group), path("pileup_focal_variants_${group}.tsv"),              emit: pileup_focal

    script:
    def g = group
    """
    set -euo pipefail

    Rscript /usr/local/bin/master_table.R "${g}"

    Rscript /usr/local/bin/binom_ggplots.R variant_master_filter_table_${g}.tsv "${g}" binom_filter_plots_${g}.pdf

    Rscript /usr/local/bin/vaf_hexbin_plots.R variant_master_filter_table_${g}.tsv "${g}" vaf_hexbin_plots_${g}.pdf

    Rscript /usr/local/bin/pileup_metric_plots.R "${g}" pileup_metric_plots_${g}.pdf

    Rscript /usr/local/bin/pileup_bppos_plots.R "${g}" variant_master_filter_table_${g}.tsv pileup_bppos_plots_${g}.pdf

    pdfunite \
        variant_filter_plot_${g}.pdf \
        binom_filter_plots_${g}.pdf \
        vaf_hexbin_plots_${g}.pdf \
        pileup_metric_plots_${g}.pdf \
        pileup_bppos_plots_${g}.pdf \
        Patient_filter_report_combined_${g}.pdf

    Rscript /usr/local/bin/vcf_annotation_table.R "${g}"

    # ── Subset res_pileup_all*.tsv files to focal (VEP-PASS) variants ─────────
    # Extract focal VariantId list (CHROM_POS_REF_ALT) from the annotation table.
    tail -n +2 vcf_annotation_table_${g}.tsv | cut -f1 | sort > focal_variant_ids_${g}.txt

    echo "Focal variants: \$(wc -l < focal_variant_ids_${g}.txt)"

    first_pileup=\$(ls res_pileup_all_group_${g}_*.tsv | head -1)
    head -1 "\${first_pileup}" > pileup_focal_variants_${g}.tsv

    # Three-lookup awk pass:
    #   Lookup 1 — ALT rows: VariantId already matches focal id; print as-is.
    #   Lookup 2 — REF rows: col6=="REF"; same CHROM+POS as a focal variant; remap VariantId
    #              and emit one row per focal variant at that position.
    #   Lookup 3 — Non-focal ALT rows: col6 is a non-REF allele but the row is at a focal
    #              CHROM+POS without matching a focal VariantId (e.g. sample carries A_G at
    #              a position where the focal variant is A_T). Zero cols 7-10 (allele-specific
    #              counts are for the wrong allele), set ALT to "REF" so the dedup step treats
    #              this as a lower-priority fallback — preserving position-level depth while
    #              letting a true ALT row win if one exists.
    #   pos_to_vid stores a pipe-delimited list so multiple focal variants sharing the same
    #   CHROM+POS all receive the real depth (fixes the single-value overwrite bug).
    awk -F'\\t' -v OFS='\\t' '
        NR==FNR {
            alt_ids[\$1] = 1
            n = split(\$1, a, "_")
            pos = a[n-2]
            chrom = a[1]; for (j=2; j<=n-3; j++) chrom = chrom "_" a[j]
            key = chrom "_" pos
            pos_to_vid[key] = (key in pos_to_vid) ? pos_to_vid[key] "|" \$1 : \$1
            next
        }
        FNR==1 { next }
        \$2 in alt_ids { print; next }
        (\$3 "_" \$4) in pos_to_vid && \$6 == "REF" {
            nvids = split(pos_to_vid[\$3 "_" \$4], vids, "|")
            for (vi = 1; vi <= nvids; vi++) {
                \$2 = vids[vi]
                print
            }
            next
        }
        (\$3 "_" \$4) in pos_to_vid {
            \$6 = "REF"; \$7 = 0; \$8 = 0; \$9 = 0; \$10 = 0
            nvids = split(pos_to_vid[\$3 "_" \$4], vids, "|")
            for (vi = 1; vi <= nvids; vi++) {
                \$2 = vids[vi]
                print
            }
        }
    ' focal_variant_ids_${g}.txt \
      res_pileup_all_group_${g}_*.tsv \
      >> pileup_focal_variants_${g}.tsv

    echo "ALT+REF pileup rows (pre-dedup): \$(tail -n +2 pileup_focal_variants_${g}.tsv | wc -l)"

    # Deduplicate: for each (sample x VariantId) pair, keep the ALT row; fall back to REF row
    # only when no ALT row exists. This prevents het sites producing two rows per pair.
    awk -F'\\t' -v OFS='\\t' '
        NR==1 { hdr=\$0; next }
        \$6 != "REF" { alt[\$1 SUBSEP \$2] = \$0; next }
                     { ref[\$1 SUBSEP \$2] = \$0 }
        END {
            print hdr
            for (k in alt) print alt[k]
            for (k in ref) if (!(k in alt)) print ref[k]
        }
    ' pileup_focal_variants_${g}.tsv > pileup_focal_dedup_${g}.tsv
    mv pileup_focal_dedup_${g}.tsv pileup_focal_variants_${g}.tsv

    echo "ALT+REF pileup rows (post-dedup): \$(tail -n +2 pileup_focal_variants_${g}.tsv | wc -l)"

    # Derive full sample list from all raw pileup files
    awk -F'\\t' 'FNR>1 {print \$1}' res_pileup_all_group_${g}_*.tsv | sort -u > all_samples_${g}.txt

    # Generate synthetic 0-rows for (sample x focal_variant) combos absent from the pileup.
    # Cols 7-22 = integer counts → 0; cols 23+ = ratios/filters/Verdict → NA.
    awk -F'\\t' -v OFS='\\t' '
        FILENAME==ARGV[1] && FNR==1 {
            ncols = NF
            filler = "0"
            for (i=8;  i<=22;    i++) filler = filler OFS "0"
            for (i=23; i<=ncols; i++) filler = filler OFS "NA"
            next
        }
        FILENAME==ARGV[2]            { samples[\$1]=1; next }
        FILENAME==ARGV[3]            { focal[\$1]=1;   next }
        FILENAME==ARGV[4] && FNR>1   { present[\$1 SUBSEP \$2]=1; next }
        END {
            for (smp in samples) {
                for (vid in focal) {
                    if (!((smp SUBSEP vid) in present)) {
                        n = split(vid, a, "_")
                        alt=a[n]; ref=a[n-1]; pos=a[n-2]
                        chrom=a[1]; for(j=2;j<=n-3;j++) chrom=chrom"_"a[j]
                        print smp, vid, chrom, pos, ref, alt, filler
                    }
                }
            }
        }
    ' pileup_focal_variants_${g}.tsv \
      all_samples_${g}.txt \
      focal_variant_ids_${g}.txt \
      pileup_focal_variants_${g}.tsv \
      >> pileup_focal_variants_${g}.tsv

    focal_rows=\$(tail -n +2 pileup_focal_variants_${g}.tsv | wc -l)
    echo "Complete pileup rows (ALT+REF+synthetic): \${focal_rows}"
    """
}
