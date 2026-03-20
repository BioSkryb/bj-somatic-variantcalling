nextflow.enable.dsl=2
params.timestamp = ""

// Merges per-sample signature cosine-similarity summaries produced by
// FILTER_VARIANTS_BY_SIG_PROBABILITY into a single cohort-level
// Samples × SBS cosine-similarity matrix.
//
// Samples without a given SBS detected as active receive NA in the matrix.
process MERGE_SIG_COSINE_SIMILARITIES {
    tag "cohort"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", mode: 'copy', enabled: "$enable_publish"

    input:
    path(filter_summaries)
    val(publish_dir)
    val(enable_publish)

    output:
    path("cohort_signature_cosine_similarities.tsv"), emit: cosine_matrix

    script:
    """
    python /scripts/merge_sig_cosine_similarities.py \\
        --input_summaries ${filter_summaries} \\
        --output          cohort_signature_cosine_similarities.tsv
    """
}
