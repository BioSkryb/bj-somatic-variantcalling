nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_VARIANT_FILTER_PROVENANCE {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/CUSTOM_VARIANT_FILTER_PROVENANCE/", enabled: "$enable_publish"
    
    input:
    tuple val(group), path(binom_beta_binom), path(dfs_gt), path(dfs_qc), path(sequoia_end)
    val(publish_dir)
    val(enable_publish)

  
    output:
    path("*")

    script:
    """

    touch done.txt;

    """
}
