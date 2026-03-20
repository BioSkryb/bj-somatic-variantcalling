nextflow.enable.dsl=2
params.timestamp = ""

// All-vs-all phylogenetic tree comparison across the 8 filtering schemes.
// Processes snv / indel / both in a single R invocation and produces:
//   - Per-type sub-directories with TSV tables and individual PDF heatmaps
//   - One master PDF combining all plots with SNV / INDEL / BOTH section headers
// Metrics: cophenetic correlation, normalised RF, branch score, Baker's Gamma,
//          and all-vs-all tanglegrams.
process TREES_COMPARE_SIMILARITIES {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), path(snv_outputs), path(indel_outputs), path(both_outputs)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group),
          path("tree_comparison_*"),
          path("${group}_tree_comparison_master.pdf"),
          emit: comparison_results

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Sentinel so the output glob always resolves even if all types are skipped
    touch tree_comparison_no_results

    Rscript /usr/local/bin/rscript_trees_compare_similarities.R \\
        --donor_id   "${group}" \\
        --tree_types "snv,indel,both" \\
        --input_dir  "." \\
        --output_dir "."
    """
}
