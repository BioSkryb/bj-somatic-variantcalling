nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_MERGE_TABLES_FILTER_SPLIT_BY_CHR {
    tag "${group}_${chr}"
    publishDir "${publish_dir}_${params.timestamp}/CUSTOM_MERGE_TABLES_FILTER_SPLIT_BY_CHR/", enabled: "$enable_publish"

    input:
    tuple val(chr), val(group), path(res_tables)
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val(group), val(chr), path("df_grouplevel_pileup_group_${group}_chr_${chr}.tsv")

    script:

    """
    cat res_* | grep -P "\\t${chr}\\t"  > df_grouplevel_pileup_group_${group}_chr_${chr}.tsv
    """
}