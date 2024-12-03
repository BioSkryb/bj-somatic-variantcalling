nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_BAM_GROUP_PILEUP_CHR {
    tag "CUSTOM_BAM_GROUP_PILEUP_CHR_${sample_name}_${chr}"
    publishDir "${publish_dir}_${params.timestamp}/CUSTOM_BAM_PILEUP_FILTER/", enabled: "$enable_publish"
    
    input:
    tuple val(group), val(chr), val(sample_name), path(input_bam), val(par_file), path(list_pos)
    path(reference)
    val(publish_dir)
    val(enable_publish)

  
    output:
    
    tuple val(sample_name), val(chr), val(par_file), path("pileup_mq0_bq0_group${group}_${sample_name}_${chr}_${par_file}.txt"), val(group)

    
    script:
    
    """
    
    samtools view -h ${input_bam[0]} ${chr} |  samtools mpileup - -f ${reference}/genome.fa --output-BP --disable-overlap-removal --count-orphans --min-BQ 0 --min-MQ 0 --output-MQ -o pileup_mq0_bq0_group${group}_${sample_name}_${chr}_${par_file}.txt -l ${list_pos}
 
  

    """
}
