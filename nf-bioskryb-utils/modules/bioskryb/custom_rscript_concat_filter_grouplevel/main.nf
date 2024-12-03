nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_RSCRIPT_CONCAT_FILTER_GROUPLEVEL {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    
    input:
    tuple val(group), path(df_files)
    path (rscript)
    val (cutoff_binomial)
    val (cutoff_beta)
    val( publish_dir )
    val( enable_publish )
  
    output:
    tuple val(group), path("res_usable_variants_${group}.tsv")

    script:
    """
    
    
    Rscript ${rscript} ${cutoff_binomial} ${cutoff_beta} 
    
    mv res.tsv res_usable_variants_${group}.tsv
    
    
    """
    
}
