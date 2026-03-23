nextflow.enable.dsl=2
params.timestamp = ""

// Adapted from CREATE_ADO_TABLE (modules/bioskryb/ado/create_ado_table/main.nf).
// The original module works with per-sample gVCFs and a separate baseline het-sites VCF.
// Here the input is already a single-sample germline VCF (output of
// SUBSET_MERGED_VCF_HIGH_CONFIDENCE_GERMLINE_FROM_STATS) which:
//   • contains only confirmed germline het sites → no additional het-site reference needed
//   • is a regular VCF (not gVCF) → gVCF conversion is not required
// Allele balance is computed directly from the FORMAT/AD and FORMAT/DP fields.
//
// Output df_ADO_* TSV format (no header, 10 columns):
//   CHROM  POS  REF  ALT  ALT  1  AD_ALT  DP  FREQ  PROVENANCE
// Columns 1–9 match what SUMMARIZE_ADO_INTERVALS expects.
// Column 10 (PROVENANCE) records the germline filter set ("stats", "vep", or "bulk")
// so downstream summaries can be split per filter type.

process CREATE_ADO_TABLE_FROM_GERMLINE_VCF {
    tag "${sample_name}"
    container 'quay.io/biocontainers/bcftools:1.14--h88f3f91_0'
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(sample_name), path(vcf), path(tbi)
    val(sample_prop)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(sample_name), path("df_ADO_${sample_name}.tsv"), emit: ado_table

    script:
    // Provenance = the trailing suffix of sample_name added by the pipeline wiring:
    // "SAMPLE_stats", "SAMPLE_vep", or "SAMPLE_bulk" → "stats" / "vep" / "bulk".
    def provenance = sample_name.tokenize('_').last()
    """
    set -euo pipefail

    # Extract allele balance at every germline het site.
    # %AD{0} = REF depth, %AD{1} = ALT depth (first ALT; germline sites are biallelic).
    # FREQ = AD_ALT / DP  (0 when DP = 0, i.e. complete allele drop-out).
    # Column 10 records the provenance label for downstream per-set summaries.
    bcftools query \\
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%AD{1}]\\t[%DP]\\n' \\
        ${vcf} \\
    | awk -v OFS='\\t' -v prop=${sample_prop} -v prov="${provenance}" \\
        'BEGIN { srand() }
        rand() < prop {
            ad_alt = \$5 + 0
            dp     = \$6 + 0
            freq   = (dp > 0) ? ad_alt / dp : 0
            print \$1, \$2, \$3, \$4, \$4, 1, ad_alt, dp, freq, prov
        }' \\
    > df_ADO_${sample_name}.tsv

    echo "[CREATE_ADO_TABLE_FROM_GERMLINE_VCF] ${sample_name}: \$(wc -l < df_ADO_${sample_name}.tsv) sites"
    """
}
