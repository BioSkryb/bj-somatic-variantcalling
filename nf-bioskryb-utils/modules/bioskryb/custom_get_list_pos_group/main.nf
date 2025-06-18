nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_GET_LIST_POS_GROUP {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/CUSTOM_GET_LIST_POS_GROUP/", enabled: "$enable_publish"
    
    input:
    tuple val(group), path(res_tables)
    val(publish_dir)
    val(enable_publish)

  
    output:
    tuple val(group), path("list_pos_variant_${group}.txt")

    
    script:
    """
    
    ls ;
  

    cat df_query_*.tsv |   cut -f1  |cut -d "_" -f1,2 | sort -u | awk -v OFS="\\t" -v FS="_" '{print \$1,\$2}' > list_pos_variant_${group}.txt


    """
}
