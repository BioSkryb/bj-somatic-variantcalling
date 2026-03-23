nextflow.enable.dsl=2
params.timestamp = ""

// Consumes all per-provenance merged_ADO_*.tsv and merged_ADO_summary_*.tsv
// files produced by CONCAT_ADO_STATS / CONCAT_ADO_VEP / CONCAT_ADO_BULK and
// generates two publication-quality plots:
//
//   ADO_germline_dist.png       — allele-frequency distribution per germline
//                                 filter set: per-sample points + red mean
//                                 line/square, faceted by provenance
//   ADO_germline_summary.png    — per-sample ADO proportion in [0.2–0.8],
//                                 dashed lines connecting the same cell
//                                 across filter sets
//   ADO_germline_comparison.png — both panels combined side-by-side
//
// The R code avoids dollar-sign column access (uses [[]]) and regex end-anchors
// so the Nextflow triple-quoted string does not require excessive escaping.

process PLOT_ADO_GERMLINE_COMPARISON {
    tag "ado_germline_plots"
    container '597246834581.dkr.ecr.us-east-1.amazonaws.com/miscellaneous:custom_snp_somatic_filter_sequoia_feb2026'
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    path(merged_ado_tables)       // merged_ADO_{stats,vep,bulk}.tsv
    path(merged_ado_summaries)    // merged_ADO_summary_{stats,vep,bulk}.tsv
    val(publish_dir)
    val(enable_publish)

    output:
    path("ADO_germline_dist.png"),        emit: dist_plot
    path("ADO_germline_summary.png"),     emit: summary_plot
    path("ADO_germline_comparison.png"),  emit: combined_plot
    path("ADO_germline_summary.tsv"),     emit: summary_table

    script:
    """
    Rscript /usr/local/bin/plot_ado_germline_comparison.R
    """
}
