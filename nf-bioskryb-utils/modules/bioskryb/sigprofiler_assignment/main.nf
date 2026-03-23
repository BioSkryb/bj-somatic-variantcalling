nextflow.enable.dsl=2
params.timestamp = ""

// Per-sample COSMIC signature assignment using SigProfilerAssignment.
// Each germline-filtered VCF is decompressed and processed individually,
// producing a per-sample results directory with signature activity scores.
process SIGPROFILER_ASSIGNMENT {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", mode: 'copy', enabled: "$enable_publish"

    input:
    tuple val(sample_name), path(vcf), path(vcf_tbi)
    val(genome_build)        // e.g. "GRCh38"
    val(context_type)        // e.g. "96"
    val(exome)                             // bool – true for exome, false for WGS
    val(export_probabilities)              // bool – export per-signature probabilities per sample
    val(export_probabilities_per_mutation) // bool – export per-mutation signature probabilities
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(sample_name), path("${sample_name}_sig_assignment"), emit: assignment_output

    script:
    def exome_flag                             = exome                             ? '--exome'                             : ''
    def export_probabilities_flag              = export_probabilities              ? '--export_probabilities'              : ''
    def export_probabilities_per_mutation_flag = export_probabilities_per_mutation ? '--export_probabilities_per_mutation' : ''
    """
    # cosmic_fit requires a directory of VCF files, not a single file path.
    mkdir -p vcf_input
    bgzip -d -c ${vcf} > vcf_input/${sample_name}.vcf

    python /scripts/run_sig_assignment.py \\
        --input_dir    vcf_input \\
        --output_dir   ${sample_name}_sig_assignment \\
        --genome_build ${genome_build} \\
        --context_type ${context_type} \\
        ${exome_flag} \\
        ${export_probabilities_flag} \\
        ${export_probabilities_per_mutation_flag}
    """
}
