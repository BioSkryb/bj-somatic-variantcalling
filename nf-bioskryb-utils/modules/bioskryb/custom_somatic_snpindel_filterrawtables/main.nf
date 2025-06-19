nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES {
    tag "${group}_${chr}"
    publishDir "${publish_dir}_${params.timestamp}/CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES/", enabled: "$enable_publish"
    
    input:
    tuple val(group),val(chr), path(res_tables), path(chosen_variants)
    val(threshold_as)
    val(threshold_clipped)
    val(threshold_prop_bp_under)  
    val(threshold_prop_bp_upper)  
    val(threshold_sd_indiv)
    val(threshold_mad_indiv)
    val(threshold_sd_both)
    val(threshold_mad_both)
    val(threshold_sd_extreme)
    val(threshold_mad_extreme)
    val(threshold_prop_cells_goodcov)
    val(threshold_good_coverage)
    val(threshold_good_variant_support)
    val(disable_qc)
    val(publish_dir)
    val(enable_publish)

  
    output:
    tuple val(group), val(chr), path("Mat_NV_${group}_${chr}.tsv"), path("Mat_NR_${group}_${chr}.tsv"), emit: tabs
    tuple val(group), val(chr), path("df_passed*"), emit:df_pass
    tuple val(group), val(chr), path("res_pileup_all_group_${group}_${chr}.tsv"), emit:pileup

    script:
    
    """

    echo -e "Concatenating tables ...";

    find . -name "res_grouplevel_pileup*" | sort -V > list_files.txt

    echo -e "SampleId\\tVariantId\\tCHROM\\tPOS\\tREF\\tALT\\tNUM_FRAGMENTS_ALLQ_MQ_BQ_F\\tNUM_FRAGMENTS_ALLQ_MQ_BQ_R\\tNUM_FRAGMENTS_HQ_MQ_BQ_F\\tNUM_FRAGMENTS_HQ_MQ_BQ_R\\tMEDIAN_AS_VARIANT_READS\\tPROP_BASES_CLIPPED\\tPROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_F\\tPROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_R\\tSD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F\\tSD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R\\tMAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F\\tMAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R\\tNUM_FRAGMENTS_ALLQ_POSITION\\tNUM_FRAGMENTS_HQ_POSITION\\tPROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_F\\tPROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_R" > res_end.tsv;

    cat list_files.txt | while read mfile; do cat \${mfile} >> res_end.tsv;done

    head -n1 res_end.tsv > df_raw_variants.tsv
    
    cat res_end.tsv | grep -v VariantId | grep -Pv "\\tREF\\t" >> df_raw_variants.tsv

    echo -e "Filtering ... ";

    Rscript /usr/local/bin/rscript_2.create_tabnr_tabnv.R ${threshold_as} ${threshold_clipped} ${threshold_prop_bp_under} ${threshold_prop_bp_upper} ${threshold_sd_indiv} ${threshold_mad_indiv} ${threshold_sd_both} ${threshold_mad_both} ${threshold_sd_extreme} ${threshold_mad_extreme} ${threshold_good_coverage} ${threshold_prop_cells_goodcov} ${threshold_good_variant_support} ${disable_qc}

    echo -e "Subsetting matrices ... ";

    head -n1 Tab_NR.tsv > Mat_NR_${group}_${chr}.tsv

    tail -n +2 Tab_NR.tsv > body.tsv

    cat ${chosen_variants} | grep "^${chr}_" > mvariants.txt

    awk -v OFS="\\t" -v FS="\\t" 'NR == FNR {  a[\$0]; next }{if(\$1 in a){ print \$0}}' mvariants.txt body.tsv >> Mat_NR_${group}_${chr}.tsv

    cat Mat_NR_${group}_${chr}.tsv | tail -n +2 | cut -f1 > vfound.txt

    awk -v OFS="\\t" -v FS="\\t" 'NR == FNR {  a[\$0]; next }{if(!(\$1 in a)){ print\$0}}' vfound.txt mvariants.txt > df_passed_PRESENTVCF_NOTINBAM_${group}_${chr}.tsv

    cat body.tsv  | cut -f1 > foundbam.txt

    awk -v OFS="\\t" -v FS="\\t" 'NR == FNR {  a[\$0]; next }{if(!(\$1 in a)){ print \$0}}' mvariants.txt foundbam.txt > df_passed_PRESENTBAM_NOTINVCF_${group}_${chr}.tsv

    head -n1 Tab_NV.tsv > Mat_NV_${group}_${chr}.tsv

    tail -n +2 Tab_NV.tsv > body.tsv

    awk -v OFS="\\t" -v FS="\\t" 'NR == FNR {  a[\$0]; next }{if(\$1 in a){ print \$0}}' ${chosen_variants} body.tsv >> Mat_NV_${group}_${chr}.tsv

    echo -e "Renaming files to output ...";
    mv df_passed_AS.tsv df_passed_AS_${group}_${chr}.tsv

    mv df_passed_propclipped.tsv df_passed_propclipped_${group}_${chr}.tsv

    mv df_passed_BPPOS.tsv df_passed_BPPOS_${group}_${chr}.tsv 

    cat Mat_NV_${group}_${chr}.tsv | tail -n +2 | cut -f1 > df_passed_DEPTH_${group}_${chr}.tsv

    mv res_end.tsv res_pileup_all_group_${group}_${chr}.tsv

    """
}
