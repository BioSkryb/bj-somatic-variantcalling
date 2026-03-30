nextflow.enable.dsl=2
params.timestamp = ""

// Compute per-sample mutational signature coverage:
//   InputVariants                    = total variants in the mutsig VCF (SNVs + INDELs)
//   NumberVariantsAssignedMutSig     = sum of all SigProfilerAssignment activity counts
//                                      (= SNVs successfully processed by SigProfiler)
//   ProportionVariantsAssignedMutSig = NumberVariantsAssignedMutSig / InputVariants
//
// Input:  {sample}_mutsig.vcf.gz from SUBSET_ANNOTATED_VCFS_FOR_MUTSIG
//         {sample}_sig_assignment/ directory from SIGPROFILER_ASSIGNMENT
// Output: mutsig_coverage_{sample}.tsv — single data row
process COMPUTE_MUTSIG_COVERAGE {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", mode: 'copy', enabled: "$enable_publish"

    input:
    tuple val(sample_name), path(mutsig_vcf), path(mutsig_tbi), path(assignment_dir)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(sample_name), path("mutsig_coverage_${sample_name}.tsv"), emit: coverage_tsv

    script:
    """
    set -euo pipefail

    # Total variants passed to SigProfiler (SNVs + INDELs)
    input_vars=\$(bcftools view -H ${mutsig_vcf} | wc -l)

    # Sum all signature activity counts (header row 1; data row 2; columns 2..NF)
    # Assignment_Solution_Activities.txt may sit in a subdirectory of the assignment dir
    acts_file=\$(find ${assignment_dir} -name "Assignment_Solution_Activities.txt" | head -1)
    assigned=\$(awk 'NR==2 {s=0; for(i=2; i<=NF; i++) s+=\$i; print s}' "\${acts_file}")

    # Proportion (guard against empty VCF)
    prop=\$(awk -v inp="\${input_vars}" -v asgn="\${assigned}" \
        'BEGIN { print (inp > 0) ? asgn/inp : 0 }')

    printf "SampleID\\tInputVariants\\tNumberVariantsAssignedMutSig\\tProportionVariantsAssignedMutSig\\n" \
        > mutsig_coverage_${sample_name}.tsv
    printf "%s\\t%s\\t%s\\t%s\\n" "${sample_name}" "\${input_vars}" "\${assigned}" "\${prop}" \
        >> mutsig_coverage_${sample_name}.tsv
    """
}
