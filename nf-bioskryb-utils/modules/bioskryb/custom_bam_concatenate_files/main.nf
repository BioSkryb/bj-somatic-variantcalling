nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_BAM_CONCATENATE_FILES {
    tag "CUSTOM_BAM_CONCATENATE_FILES"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/alignment/subsample", enabled:"$enable_publish"


    input:
    tuple val(output_name), path(bam_files)
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val("${output_name}"), path("${output_name}.bam"), path("${output_name}.bam.bai"), emit: bam

    script:

    """
    
    samtools merge --threads $task.cpus temp_merged.bam `find . -name "*.bam"`
    
    samtools sort --threads $task.cpus -O BAM -o ${output_name}.bam temp_merged.bam
    samtools index ${output_name}.bam
    
    """
}