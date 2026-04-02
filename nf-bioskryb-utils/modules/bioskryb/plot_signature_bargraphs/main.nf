nextflow.enable.dsl=2
params.timestamp = ""

// Generates signature activity bar graphs from the zero-filtered and
// cosine-filtered activity matrices produced by FILTER_ACTIVITIES_BY_COSINE,
// combined with the cosine similarity matrix from MERGE_SIG_COSINE_SIMILARITIES.
//
// Outputs:
//   signature_bargraphs_combined.png  — static three-panel PNG (top: proportion assigned
//                                       mutational signature, Panel A: unfiltered,
//                                       Panel B: cosine filtered), shared x-axis sample order
//                                       from Bray-Curtis clustering.
//   signature_bargraphs_interactive.html — D3-based interactive version with
//                                          hover tooltips (cosine sim per segment)
//                                          and click-to-highlight legend.
process PLOT_SIGNATURE_BARGRAPHS {
    tag "cohort"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", mode: 'copy', enabled: "$enable_publish"

    input:
    path(zero_filtered_activities)    // zero_filtered_signature_activities.tsv
    path(cosine_filtered_activities)  // cosine_filtered_signature_activities.tsv
    path(cosine_similarities)         // cohort_signature_cosine_similarities.tsv
    path(mutsig_coverage)             // mutsig_coverage.tsv
    val(publish_dir)
    val(enable_publish)

    output:
    path("signature_bargraphs_combined.png"),     emit: bargraph_png
    path("signature_bargraphs_interactive.html"), emit: bargraph_html

    script:
    """
    Rscript /usr/local/bin/plot_signature_bargraphs.R \\
        ${zero_filtered_activities} \\
        ${cosine_filtered_activities} \\
        ${mutsig_coverage}

    Rscript /usr/local/bin/generate_interactive_bargraph.R \\
        ${zero_filtered_activities} \\
        ${cosine_filtered_activities} \\
        ${cosine_similarities}
    """
}
