nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_RSCRIPT_FILTER_GROUPLEVEL_PILEUP_SAMPLE {
    tag "${group}_${chr}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    
    input:
    tuple val(group), val(gender), val(chr), path(df_files)
    path (rscript)
    val (cutoff_depth)
    val (cutoff_nsupvar)
    val (cutoff_vaf)
    val (cutoff_filter_noise_prev)
    val (cutoff_presence_prev)
    val( publish_dir )
    val( enable_publish )
  
    output:
    tuple val(group), val(chr), path("res_grouplevel_filtered_pileup_${chr}_${group}.tsv")

    script:
    """
    
    Rscript ${rscript} ${cutoff_depth} ${cutoff_nsupvar} ${cutoff_vaf} ${cutoff_filter_noise_prev} ${cutoff_presence_prev} ${gender}
    
    mv res.tsv res_grouplevel_filtered_pileup_${chr}_${group}.tsv
    
    
    """
    
}
