nextflow.enable.dsl=2
params.timestamp = ""

// Merge per-sample mutsig_coverage_*.tsv files produced by COMPUTE_MUTSIG_COVERAGE
// into a single cohort-level TSV with one row per sample.
// Output columns: SampleID | InputVariants | NumberVariantsAssignedMutSig | ProportionVariantsAssignedMutSig
process MERGE_MUTSIG_COVERAGE {
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", mode: 'copy', enabled: "$enable_publish"

    input:
    path(coverage_files)   // collected list of per-sample mutsig_coverage_*.tsv
    val(publish_dir)
    val(enable_publish)

    output:
    path("mutsig_coverage.tsv"), emit: merged_coverage

    script:
    """
    set -euo pipefail
    head -1 \$(ls mutsig_coverage_*.tsv | head -1) > mutsig_coverage.tsv
    for f in mutsig_coverage_*.tsv; do tail -n +2 "\$f"; done >> mutsig_coverage.tsv
    """
}
