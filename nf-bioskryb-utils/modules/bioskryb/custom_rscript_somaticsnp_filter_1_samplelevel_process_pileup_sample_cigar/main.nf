nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_RSCRIPT_SOMATICSNP_FILTER_1_SAMPLELEVEL_PROCESS_PILEUP_SAMPLE_CIGAR {
    tag "${sample_name}_${chr}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"


    input:
    tuple val(group), val(chr), val(sample_name), path(pileup_cigars_file), path(chosen_variants)
    val(threshold_mq)
    val(threshold_bq)
    val(threshold_bp)
    val(num_lines_read_pileup)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group),val(chr), path("res_grouplevel_pileup_group_${group}_sample_${sample_name}_chr_${chr}.tsv")

    script:
    """
    echo -e "Subsetting variants ...";
    cat ${chosen_variants} | grep "^${chr}_"  > cvariants.txt
    cat cvariants.txt | cut -d "_" -f2 | sort -uV > cpositions.txt
    
    awk -v OFS="\\t" 'NR == FNR {  a[\$1]=1; next }{if( \$2 in a ){print \$0}}' cpositions.txt ${pileup_cigars_file} > pileup_subset.tsv
    echo -e "Launching QC script ...";
    Rscript /usr/local/bin/rscript_1.samplelevel_process_pileup_sample_cigars.R pileup_subset.tsv ${sample_name} ${threshold_mq} ${threshold_bq} ${threshold_bp} ${num_lines_read_pileup}
    mv res.tsv res_grouplevel_pileup_group_${group}_sample_${sample_name}_chr_${chr}.tsv
    
    """

}