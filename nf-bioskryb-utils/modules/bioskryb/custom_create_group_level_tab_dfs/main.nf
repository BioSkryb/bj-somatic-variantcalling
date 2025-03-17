nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_CREATE_GROUP_LEVEL_TAB_DFS {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/CUSTOM_CREATE_GROUP_LEVEL_TAB_DFS/", enabled: "$enable_publish"

    input:
    tuple val(group), path(res_tables_nv),  path(res_tables_nr)
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val(group), path("Mat_NV_${group}.tsv"), path("Mat_NR_${group}.tsv"), emit: tabs

    script:

    """
    ls Mat_NR* | sort -V > list_nr.txt
    ls Mat_NV* | sort -V > list_nv.txt
    ls Mat_N* | while read file; do head -n1 \${file} | tr "\\t" "\\n" | tail -n+2;done | sort -uV > all_samples.txt
    Rscript /usr/local/bin/rscript_3.create_tabnr_tabnv_group.R
    mv Tab_NR.tsv Mat_NR_${group}.tsv
    mv Tab_NV.tsv Mat_NV_${group}.tsv
    """
}