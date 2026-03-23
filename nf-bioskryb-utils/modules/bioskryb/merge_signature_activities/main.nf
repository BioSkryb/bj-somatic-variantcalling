nextflow.enable.dsl=2
params.timestamp = ""

// Merges per-sample Assignment_Solution_Activities.txt files produced by
// SIGPROFILER_ASSIGNMENT into a single cohort-level SamplesXSignatures matrix.
// All assignment output directories are staged into the work directory and
// processed by merge_signature_activities.py, which unions signature columns
// and fills absent values with 0.
process MERGE_SIGNATURE_ACTIVITIES {
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", mode: 'copy', enabled: "$enable_publish"

    input:
    path(assignment_dirs)   // collected list of per-sample SigProfilerAssignment output directories
    val(publish_dir)
    val(enable_publish)

    output:
    path("merged_signature_activities.txt"), emit: merged_activities

    script:
    """
    python /scripts/merge_signature_activities.py \\
        --input_dirs ${assignment_dirs} \\
        --output     merged_signature_activities.txt
    """
}
