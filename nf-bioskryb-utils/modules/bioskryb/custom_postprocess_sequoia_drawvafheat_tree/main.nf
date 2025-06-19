nextflow.enable.dsl=2
params.timestamp = ""

process POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    

    
    input:
    tuple val(group), path(file_placed), path(file_nv), path(file_nr), path(file_tree)
    path (df_gt_files)
    val( publish_dir )
    val( enable_publish )
  
    output:
    tuple val(group), path("res_composition_tree_vafheatmap_${group}.pdf"), path("res_composition_tree_digitalheatmap_${group}.pdf"), path("res_figures_${group}.RDS")

    script:
    """
    
    ls df_gt* | while read file; 
    do
        name=`echo \${file} | sed 's/df_gt_//' | sed 's/\\.tsv//'`;
        cat \${file} | sed "s|^|\${name}\\t|";
    done > df_all_gt.tsv

    cat ${file_nr} | cut -d " " -f1 | tail -n+2  | sed 's/\\"//g' > ids.txt

    awk -v OFS="\\t" -v FS="\\t" 'NR == FNR {  a[\$0]; next }{if(\$2 in a){ print \$0}}' ids.txt df_all_gt.tsv  > df_all_gt_chosen.tsv

    Rscript /usr/local/bin/rscript_4.postprocess_vaf_drawtree_heatmap_intnodes.R ${file_placed} ${file_nv} ${file_nr} ${file_tree} df_all_gt_chosen.tsv
    
    mv res_composition.pdf res_composition_tree_vafheatmap_${group}.pdf

    mv res_composition_digital.pdf res_composition_tree_digitalheatmap_${group}.pdf

    mv res_figures.RDS res_figures_${group}.RDS
    
    """
    
}
