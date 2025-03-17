nextflow.enable.dsl=2
params.timestamp = ""

process CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"

    input:
    tuple val(group), path(files)
    val(cutoff_binomial)
    val(cutoff_beta)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), path("res_${group}_binomial_betabinomial.tsv"), emit: res_df
    tuple val(group), path("chosen_variants_${group}.txt"), emit: chosen_variants

    script:
    """
    
    echo -e "Concatenating tables ...";
    cat filtered*.txt | grep "Depth_filter" | sed "s|^|VariantId\\t|" | sort -u > res.tsv
    ls filtered*.txt | sort -V | while read file; do cat \${file} | tail -n+2 >> res.tsv  ;done
    echo -e "Launching filtering with binomial cutoff: ${cutoff_binomial} and beta-binomial cutoff: ${cutoff_beta}";
    Rscript /usr/local/bin/rscript_0.filter_binom_betabinom_tables.R res.tsv ${cutoff_binomial} ${cutoff_beta}
    mv df_verdict.txt res_${group}_binomial_betabinomial.tsv
    mv chosen_variants.txt chosen_variants_${group}.txt
    
    """

}