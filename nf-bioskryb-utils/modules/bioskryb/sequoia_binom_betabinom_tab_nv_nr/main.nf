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

    head -n1 ${tab_nv} > header_nv.tsv;

    head -n1 ${tab_nr} > header_nr.tsv;

    cat ${tab_nv} | tail -n +2 | grep -v "NON_REF" | cat header_nv.tsv - > mat_nv.tsv

    cat ${tab_nr} | tail -n +2 | grep -v "NON_REF" | cat header_nr.tsv - > mat_nr.tsv
 
    wc -l mat_nv.tsv;
    
    wc -l mat_nr.tsv;

    echo -e "Launching filtering pipeline ...";

    Rscript /usr/local/bin/sequoia_filter_variants_nvnr.R -v mat_nv.tsv -r mat_nr.tsv -n 1 --min_cov ${par_min_cov} --max_cov ${par_max_cov} --gender_passed ${gender}

    ls *.txt | while read file; do mv \${file} filtered_${group}_${chr}_\${file};done

    echo -e "Done";
    
    """
    
}
