nextflow.enable.dsl=2
params.timestamp = ""

// Extract the sorted list of sample names from a group-level VCF.
// Emits a single sample_list.txt per group; used to fan-out per-sample tasks downstream.
process LIST_SAMPLES_FROM_GROUP_VCF {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), path(vcf), path(tbi)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), path("sample_list.txt"), emit: sample_list

    script:
    """
    bcftools query -l ${vcf} | sort > sample_list.txt
    echo "Samples found: \$(wc -l < sample_list.txt)"
    """
}
