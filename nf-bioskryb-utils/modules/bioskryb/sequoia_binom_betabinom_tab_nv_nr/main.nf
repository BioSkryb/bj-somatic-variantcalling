nextflow.enable.dsl=2
params.timestamp = ""

process SEQUOIA_BINOM_BETABINOM_TAB_NV_NR {
    tag "${group}_${chr}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"

    input:
    tuple val(group), val(chr), path(tab_nv), path(tab_nr)
    val(par_min_cov)
    val(par_max_cov)
    val(gender)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group),val(chr), path("*_filtering_all.txt"), emit: df_filter

    script:
    """
    
    wc -l ${tab_nv};
    wc -l ${tab_nr};
    echo -e "Launching filtering pipeline ...";
    Rscript /usr/local/bin/sequoia_filter_variants_nvnr.R -v ${tab_nv} -r ${tab_nr} -n 1 --min_cov ${par_min_cov} --max_cov ${par_max_cov} --gender_passed ${gender}
    ls *.txt | while read file; do mv \${file} filtered_${group}_${chr}_\${file};done
    echo -e "Done";
    
    """

}