nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_RSCRIPT_SOMATICSNP_FILTER_2_CELL_LEVEL_CONCAT_FILTER_TABLE {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    
    input:
    tuple val( sample_name ), val(group), path(res_tables)
    val(cutoff_as)
    val(cutoff_prop_clipped_reads)
    val(cutoff_num_hq_frag)
    val(filter_bp_pos)
    val(cutoff_prop_bp)
    val(cutoff_sd_indiv)
    val(cutoff_mad_indiv)
    val(cutoff_sd_both)
    val(cutoff_mad_both)
    val(cutoff_sd_extreme)
    val(cutoff_mad_extreme)

    val( publish_dir )
    val( enable_publish )
  
    output:
    tuple val(group), val(sample_name), path("res_table_raw_somaticsnp_${sample_name}.tsv"), emit: res_raw
    tuple val(group), val(sample_name), path("res_table_filtered_somaticsnp_${sample_name}.tsv"), emit:res_filter
    tuple val(group), val(sample_name), path("res_*.tsv"), emit: all_tables

    script:
    """


    Rscript /usr/local/bin/rscript_2.concat_filter_celllevel_table_baminfovariants.R ${cutoff_as} ${cutoff_prop_clipped_reads} ${cutoff_num_hq_frag} ${filter_bp_pos} ${cutoff_prop_bp} ${cutoff_sd_indiv} ${cutoff_mad_indiv} ${cutoff_sd_both} ${cutoff_mad_both} ${cutoff_sd_extreme} ${cutoff_mad_extreme}
    
    mv res_raw.tsv res_table_raw_somaticsnp_${sample_name}.tsv
    
    mv res_1as_2propclip_3numhq_4posvar.tsv res_table_filtered_somaticsnp_${sample_name}.tsv
    
    """
    
}
