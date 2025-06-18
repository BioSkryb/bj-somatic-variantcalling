nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_RSCRIPT_SOMATICSNP_FILTER_4_CREATE_TABNR_TABNV {
    tag "${group}_${chr}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    

    
    input:
    tuple val(group), val(chr), path(df_files)
    val(cutoff_depth_manual)
    val(cutoff_numreads_variant_manual)
    val(cutoff_prev_na_manual)
    val(cutoff_prev_var_manual)
    val( publish_dir )
    val( enable_publish )
  
    output:
    tuple val(group), path("Tab_NR_somatic_snps_${group}_${chr}.tsv"), path("Tab_NV_somatic_snps_${group}_${chr}.tsv"), emit:raw
    tuple val(group), path("Tab_*filtered_${group}_${chr}.tsv"), emit:clean

    script:
    """
    
    Rscript /usr/local/bin/rscript_4.create_tabnr_tabns.R ${cutoff_depth_manual} ${cutoff_numreads_variant_manual} ${cutoff_prev_na_manual} ${cutoff_prev_var_manual} ${chr}
    
    
    mv Tab_NR_somatic_snps.tsv Tab_NR_somatic_snps_${group}_${chr}.tsv
    
    
    mv Tab_NV_somatic_snps.tsv Tab_NV_somatic_snps_${group}_${chr}.tsv
    
    
    
    mv Tab_NR_somatic_snps_filtered.tsv Tab_NR_somatic_snps_filtered_${group}_${chr}.tsv
    
    
    mv Tab_NV_somatic_snps_filtered.tsv Tab_NV_somatic_snps_filtered_${group}_${chr}.tsv
        
    """
    
}
