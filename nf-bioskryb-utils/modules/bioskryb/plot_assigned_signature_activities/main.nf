nextflow.enable.dsl=2
params.timestamp = ""

// Generates a stacked bar chart from the merged SamplesXSignatures activity
// matrix produced by MERGE_SIGNATURE_ACTIVITIES.
// X-axis: samples (90° tick labels); Y-axis: number of mutations.
// Bars are stacked by COSMIC signature ID and coloured with a distinct palette.
// Only signatures with at least one non-zero count are shown.
process PLOT_ASSIGNED_SIGNATURE_ACTIVITIES {
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", mode: 'copy', enabled: "$enable_publish"

    input:
    path(merged_activities)   // activities matrix TSV
    val(label)                // filename prefix, e.g. "merged", "zero_filtered", "cosine_filtered"
    val(publish_dir)
    val(enable_publish)

    output:
    path("${label}_signature_activities.png"), emit: signature_plot

    script:
    """
    python /scripts/plot_signature_activities.py \\
        --input  ${merged_activities} \\
        --output ${label}_signature_activities.png
    """
}
