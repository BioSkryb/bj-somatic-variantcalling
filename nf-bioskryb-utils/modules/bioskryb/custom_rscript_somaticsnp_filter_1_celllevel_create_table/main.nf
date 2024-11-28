nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_RSCRIPT_SOMATICSNP_FILTER_1_CELL_LEVEL_CREATE_TABLE {
    tag "${sample_name}_${chr}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    
    input:
    tuple val( sample_name ), val(chr), val(paral_id), path(pileup_file), path( df_cig ), val(group)
    val( publish_dir )
    val( enable_publish )
  
    output:
    tuple val(sample_name), val(group), path("res_table_filter_${sample_name}_${chr}_${paral_id}.tsv")

    script:
    """
    
    Rscript /usr/local/bin/rscript_1.create_celllevel_table_baminfovariants.R  1 ${pileup_file}
    
    
    cat res_flags.tsv | cut -f1 | grep -v NA_NA | sort -u > res_flags_unique.txt
    

    awk -v FS="\\t" 'NR == FNR {  a[\$0]; next }{var=\$1"_"\$2; if(var in a){ print \$0 }}' res_flags_unique.txt ${df_cig}  > df_found



    Rscript /usr/local/bin/rscript_1.create_celllevel_table_baminfovariants.R 2


    
    mv res_end.tsv res_table_filter_${sample_name}_${chr}_${paral_id}.tsv 

    
    """
    
}
