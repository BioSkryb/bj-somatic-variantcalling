nextflow.enable.dsl=2
params.timestamp = ""

process SUMMARIZE_ADO_INTERVALS {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    
    input:
    tuple val( sample_name ), path( ado_file )
    val( cov_cutoff )
    val( publish_dir )
    val( enable_publish )
  
    output:
    path("*.png"), emit: plot
    path("res_ADO_*"), emit: df_sum

    script:
    """
    
    echo -e "Processing ${ado_file} with coverage cutoff >= ${cov_cutoff}";
    Rscript /usr/local/bin/summarize_ado_intervals.R ${ado_file} ${cov_cutoff}
    """
    
}
