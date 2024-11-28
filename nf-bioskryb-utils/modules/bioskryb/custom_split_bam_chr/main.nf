nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_SPLIT_BAM_CHR {
    tag "CUSTOM_SPLIT_BAM_CHR_${sample_name}_${chr}"
    publishDir "${publish_dir}_${params.timestamp}/CUSTOM_SPLIT_BAM_CHR/", enabled: "$enable_publish"
    
    input:
    tuple val(sample_name), path(input_bam), val(group), val(chr)
    val(publish_dir)
    val(enable_publish)

  
    output:
    tuple val(sample_name), val(chr), path("${sample_name}_${chr}.bam*"), val(group), emit:bam
    tuple val(sample_name), val(chr), path("df_name_flag_cigar_${sample_name}_${chr}.tsv"), val(group), emit:df

    
    script:
    
    """
    samtools view --threads $task.cpus -b -o ${sample_name}_${chr}_unsorted.bam ${input_bam[0]} ${chr} 
    
    samtools sort --threads $task.cpus -o ${sample_name}_${chr}.bam ${sample_name}_${chr}_unsorted.bam
    
    samtools index --threads $task.cpus ${sample_name}_${chr}.bam


    samtools view --threads $task.cpus ${sample_name}_${chr}.bam | cut -f1,2,6 > df_name_flag_cigar_${sample_name}_${chr}.tsv
    

 

    """
}