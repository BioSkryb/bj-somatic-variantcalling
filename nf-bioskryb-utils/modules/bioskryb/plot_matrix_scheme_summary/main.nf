nextflow.enable.dsl=2
params.timestamp = ""

// Render a PDF table showing SNV / INDEL / Total counts for each of the
// NR/NV matrix filtering schemes produced by CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF,
// plus boxplots of per-sample variant counts at each upstream filter stage
// (pre_bulk, post_bulk, post_binom, post_vep) from CUSTOM_VARIANT_FILTER_PROVENANCE.
//
// Inputs:
//   scheme_tsv           — matrix_scheme_summary_${group}.tsv
//                            columns: scheme | NumberOfSNVs | NumberOfIndels
//   per_sample_tsv       — matrix_per_sample_summary_${group}.tsv
//                            columns: sample | scheme | NumberOfSNVs | NumberOfIndels
//   upstream_per_sample_tsv — upstream_filter_per_sample_${group}.tsv
//                            same columns; stages: pre_bulk, post_bulk, post_binom, post_vep
// Output: matrix_scheme_summary_${group}.pdf
process PLOT_MATRIX_SCHEME_SUMMARY {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), path(scheme_tsv), path(per_sample_tsv), path(upstream_per_sample_tsv)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), path("matrix_scheme_summary_${group}.pdf"), emit: scheme_summary_pdf

    script:
    """

    set -euo pipefail
    Rscript /usr/local/bin/plot_matrix_scheme_summary.R \
        "${group}" "${scheme_tsv}" "${per_sample_tsv}" "${upstream_per_sample_tsv}"

    """
}
