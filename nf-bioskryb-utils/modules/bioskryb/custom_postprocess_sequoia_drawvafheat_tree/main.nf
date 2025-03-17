nextflow.enable.dsl=2
params.timestamp = ""

process POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"



    input:
    tuple val(group), path(file_placed), path(file_nv), path(file_nr), path(file_tree)
    val( publish_dir )
    val( enable_publish )

    output:
    tuple val(group), path("res_*")

    script:
    """
    
    Rscript /usr/local/bin/rscript_4.postprocess_vaf_drawtree_heatmap_intnodes.R ${file_placed} ${file_nv} ${file_nr} ${file_tree}
    
    mv res_composition.pdf res_composition_tree_vafheatmap_${group}.pdf
    mv res_vaf.tsv res_vaf_${group}.tsv
    mv res_vaf_ptree.RDS res_vaf_ptree_${group}.RDS
    
    """

}