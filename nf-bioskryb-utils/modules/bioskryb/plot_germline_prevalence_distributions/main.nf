nextflow.enable.dsl=2
params.timestamp = ""

// Reads the long-format prevalence TSV produced by EXTRACT_GERMLINE_PREVALENCE_TABLE
// and generates a publication-quality histogram figure:
//
//   germline_prevalence_distributions_{group}.png
//       — one facet per filter (GERMLINE_FROM_STATS, VEP_AF_filter, Bulk_Fail)
//       — bars coloured by pass/fail relative to germline_prev_pct threshold
//       — vertical dashed line at the threshold
//       — variant counts and median annotated per panel

process PLOT_GERMLINE_PREVALENCE_DISTRIBUTIONS {
    tag "${group}"
    container '597246834581.dkr.ecr.us-east-1.amazonaws.com/miscellaneous:custom_snp_somatic_filter_sequoia_feb2026'
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), path(prevalence_table)
    val(germline_prev_pct)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), path("germline_prevalence_distributions_${group}.png"), emit: prevalence_plot

    script:
    """
    echo "${germline_prev_pct}" > threshold.txt

    Rscript /usr/local/bin/plot_germline_prevalence_distributions.R

    # Rename placeholder to group-specific filename
    mv germline_prevalence_distributions_PLACEHOLDER.png \
       germline_prevalence_distributions_${group}.png
    """
}
