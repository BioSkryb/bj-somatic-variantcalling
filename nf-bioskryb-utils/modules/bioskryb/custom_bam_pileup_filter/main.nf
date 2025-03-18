nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_BAM_PILEUP_FILTER {
    tag "CUSTOM_BAM_PILEUP_FILTER_${sample_name}_${chr}"
    publishDir "${publish_dir}_${params.timestamp}/CUSTOM_BAM_PILEUP_FILTER/", enabled: "$enable_publish"

    input:
    tuple val(sample_name), val(chr), val(paral_id), path(list_pos), path(input_bam), val(group)
    path(reference)
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val(sample_name), val(chr), val(paral_id), path("pileup_mq0_bq0_${sample_name}_${chr}_${paral_id}.txt")


    script:

    """
 
    samtools mpileup ${input_bam[0]}  -f ${reference}/genome.fa --output-BP --output-QNAME  --disable-overlap-removal --count-orphans --min-BQ 0 --min-MQ 0 --output-MQ --output-extra FLAG,AS -o pileup_mq0_bq0_${sample_name}_${chr}_${paral_id}.txt -l ${list_pos}
 
    """
}