nextflow.enable.dsl=2
params.timestamp = ""

// Zeros out (sample, SBS) activity counts where the per-signature cosine
// similarity is below min_cosine_threshold, indicating that the mutations
// attributed to that signature do not resemble its expected COSMIC profile.
//
// Inputs:
//   cosine_matrix – cohort_signature_cosine_similarities.tsv (Samples × SBS)
//   merged_activities – merged_signature_activities.txt (Samples × SBS counts)
//
// Outputs:
//   cosine_filtered_signature_activities.tsv – filtered count matrix
//   cosine_filter_report.tsv                 – log of zeroed (sample, SBS) pairs
process FILTER_ACTIVITIES_BY_COSINE {
    tag "cohort"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", mode: 'copy', enabled: "$enable_publish"

    input:
    path(cosine_matrix)
    path(merged_activities)
    val(min_cosine)
    val(publish_dir)
    val(enable_publish)

    output:
    path("zero_filtered_signature_activities.tsv"),   emit: zero_filtered_activities
    path("cosine_filtered_signature_activities.tsv"), emit: filtered_activities
    path("cosine_filter_report.tsv"),                 emit: filter_report
    path("cosine_filter_stats.txt"),                  emit: filter_stats

    script:
    """
    python /scripts/filter_signature_activities_by_cosine.py \\
        --cosine_matrix        ${cosine_matrix} \\
        --activities           ${merged_activities} \\
        --min_cosine           ${min_cosine} \\
        --output_zero_filtered zero_filtered_signature_activities.tsv \\
        --output               cosine_filtered_signature_activities.tsv \\
        --output_report        cosine_filter_report.tsv \\
        --output_stats         cosine_filter_stats.txt
    """
}
