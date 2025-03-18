nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_BAM_GROUP_PILEUP_CHR {
    tag "${sample_name}_${chr}_${par_file}"
    publishDir "${publish_dir}_${params.timestamp}/CUSTOM_BAM_PILEUP_FILTER/", enabled: "$enable_publish"

    input:
    tuple val(group), val(chr), val(sample_name), path(input_bam), path(df_flags), val(par_file), path(list_pos)
    path(reference)
    val(publish_dir)
    val(enable_publish)


    output:  
    tuple val(sample_name), val(chr), val(par_file), path("pileup_mq0_bq0_group${group}_${sample_name}_${chr}_${par_file}.txt"), path("df_name_flag_cigar_${group}_${sample_name}_${chr}_${par_file}.txt"), val(group)


    script:

    """
    
    samtools view --threads $task.cpus -h ${input_bam[0]} ${chr} |  samtools mpileup - -f ${reference}/genome.fa --output-BP-5 --output-QNAME --disable-overlap-removal --count-orphans --min-BQ 0 --min-MQ 0 --output-MQ --output-extra FLAG,AS -o pileup_mq0_bq0_group${group}_${sample_name}_${chr}_${par_file}.txt -l ${list_pos}
    
    cut -f8 pileup_mq0_bq0_group${group}_${sample_name}_${chr}_${par_file}.txt | tr "," "\\n" > reads_touse_a.txt
    cut -f9 pileup_mq0_bq0_group${group}_${sample_name}_${chr}_${par_file}.txt | tr "," "\\n" > reads_touse_b.txt
    paste -d "_" reads_touse_a.txt reads_touse_b.txt  | sort -u > reads_touse.txt
    awk -v OFS="\\t" -v FS="\\t" 'NR == FNR {  a[\$0]; next }{uid=\$1"_"\$2;if(uid in a){ print \$1"_"\$2,\$3 }}' reads_touse.txt ${df_flags}  > df_name_flag_cigar_${group}_${sample_name}_${chr}_${par_file}.txt
    """
}