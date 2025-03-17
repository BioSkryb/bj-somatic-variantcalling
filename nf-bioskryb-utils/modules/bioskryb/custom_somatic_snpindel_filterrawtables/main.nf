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
    val(threshold_sd_indiv)
    val(threshold_mad_indiv)
    val(threshold_sd_both)
    val(threshold_mad_both)
    val(threshold_sd_extreme)
    val(threshold_mad_extreme)
    val(threshold_prop_cells_goodcov)
    val(threshold_good_coverage)
    val(threshold_good_variant_support)
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
    echo -e "SampleId\\tVariantId\\tCHROM\\tPOS\\tREF\\tALT\\tNUM_FRAGMENTS_ALLQ_MQ_BQ_F\\tNUM_FRAGMENTS_ALLQ_MQ_BQ_R\\tNUM_FRAGMENTS_HQ_MQ_BQ_F\\tNUM_FRAGMENTS_HQ_MQ_BQ_R\\tMEDIAN_AS_VARIANT_READS\\tPROP_BASES_CLIPPED\\tPROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_F\\tPROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_R\\tSD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F\\tSD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R\\tMAD_BPSTART_RAGMENTS_HQ_MQ_BQ_F\\tMAD_BPSTART_RAGMENTS_HQ_MQ_BQ_R\\tNUM_FRAGMENTS_ALLQ_POSITION\\tNUM_FRAGMENTS_HQ_POSITION" > res_end.tsv;
    cat list_files.txt | while read mfile; do cat \${mfile} >> res_end.tsv;done
    ##############################
    ### Sample-level filtering ###
    ##############################
    cat res_end.tsv | grep -v VariantId | grep -Pv "\\tREF\\t" > df_raw_variants.tsv
    #Filter variants by AS
    echo -e "Sample-level filtering: AS ... ";
    cat df_raw_variants.tsv  | awk -v OFS="\\t" -v t_as="${threshold_as}" '{if(\$11 >= t_as){print \$1,\$2}}' | sed "s|^|PASSED_AS\\t|" > df_passed_AS.tsv
    #Filter variants by proportion of clipped reads supporting the variant
    echo -e "Sample-level filtering: PropClipped ... ";
    cat df_raw_variants.tsv  | awk -v OFS="\\t" -v t_pc="${threshold_clipped}" '{if(\$12 < t_pc){print \$1,\$2}}'  | sed "s|^|PASSED_PROPCLIPPED\\t|" > df_passed_propclipped.tsv
    #Filter variants based of the variants positions in the reads
    echo -e "Sample-level filtering: VariantPosition ... ";
    cat df_raw_variants.tsv | awk -v OFS="\\t" -v t_prop="${threshold_prop_bp_under}" -v t_sd="${threshold_sd_indiv}" -v t_mad="${threshold_mad_indiv}" '{if(\$9 <2 && \$10 >1){if( (\$14 < t_prop) || ( \$16 > t_sd && \$18 > t_mad) ){print \$1,\$2}}}' > df_a.tsv
    cat df_raw_variants.tsv | awk -v OFS="\\t" -v t_prop="${threshold_prop_bp_under}" -v t_sd="${threshold_sd_indiv}" -v t_mad="${threshold_mad_indiv}" '{if(\$9 >1 && \$10 <2){if( (\$13 < t_prop) || ( \$15 > t_sd && \$17 > t_mad) ){print \$1,\$2}}}' > df_b.tsv
    cat df_raw_variants.tsv | awk -v OFS="\\t" '{if(\$9 >1 && \$10 >1){print \$0}}' > _temp_c.tsv
    cat _temp_c.tsv |  awk -v OFS="\\t" -v t_prop="${threshold_prop_bp_under}" -v t_sd_b="${threshold_sd_both}" -v t_mad_b="${threshold_mad_both}" -v t_sd_e="${threshold_sd_extreme}" -v t_mad_e="${threshold_mad_extreme}"  '{if( (\$13 < t_prop)  && ( ( \$15 > t_sd_b && \$17 > t_mad_b) || (\$16 > t_sd_e && \$18 > t_mad_e) ) ){print \$1,\$2}}' > df_c_1.tsv
    cat _temp_c.tsv |  awk -v OFS="\\t" -v t_prop="${threshold_prop_bp_under}" -v t_sd_b="${threshold_sd_both}" -v t_mad_b="${threshold_mad_both}" -v t_sd_e="${threshold_sd_extreme}" -v t_mad_e="${threshold_mad_extreme}" '{if( (\$14 < t_prop)  && ( ( \$16 > t_sd_b  && \$18 > t_mad_b) || (\$15 > t_sd_e && \$17 > t_mad_e) ) ){print \$1,\$2}}' > df_c_2.tsv
    #Summarize the union of the variants passin the filter of the bp position in the reads
    cat df_a.tsv df_b.tsv df_c_1.tsv df_c_2.tsv | sort -u | sed "s|^|PASSED_BPPOSINREADS\\t|" > df_passed_BPPOS.tsv
    #Variants that pass to the group level are the ones that pass the three filters in at least one sample
    cat df_passed_AS.tsv df_passed_propclipped.tsv df_passed_BPPOS.tsv | awk '{print \$2"|"\$3}' | sort | uniq -c  | awk '{if (\$1 ==3){print \$2}}' | cut -d "|" -f2 | sort -u > variants_to_use_group.txt
    ##############################
    ### Group level filtering ####
    ##############################
    echo -e "Group-level filtering: Preparation ... ";
    cat res_end.tsv | grep -v VariantId  > df_raw_variants.tsv
    #Get all the alt ids
    awk -v OFS="\\t" -v FS="\\t" 'NR == FNR {  a[\$0]; next }{if(\$2 in a){ print \$0}}' variants_to_use_group.txt df_raw_variants.tsv > df_raw_variants_group_alt.tsv
    
    #Get all the ref ids too
    cat variants_to_use_group.txt | cut -d "_" -f1,2 | sort -u > pos_to_use_group.txt
    
    cat df_raw_variants.tsv  | grep -P "\\tREF\\t" | awk -v OFS="\\t" '{print \$3"_"\$4,\$0}' > df_raw_variants_ref.tsv
    awk -v OFS="\\t" -v FS="\\t" 'NR == FNR {  a[\$0]; next }{if(\$1 in a){ print \$0}}' pos_to_use_group.txt df_raw_variants_ref.tsv | cut -f2- > df_raw_variants_group_ref.tsv
    cat df_raw_variants_group_ref.tsv df_raw_variants_group_alt.tsv > df_raw_variants_group.tsv
    #Filter based on total depth in the position. The idea is that at least x% of the cells should have coverage of at least y in the position
    num_samples=`cat df_raw_variants.tsv | cut -f1 | sort -u  | wc -l`;
    #Calculate proportion of samples required to pass
    threshold_pass_group=`echo "${threshold_prop_cells_goodcov} * \${num_samples}" | bc -l`;
    echo -e "Group-level filtering: Depth ... ";
    cut -f1,2,20 df_raw_variants_group.tsv | sort -u | awk -v t_gc="${threshold_good_coverage}" '{if(\$3 >= t_gc){print \$0}}' | cut -f2 | sort | uniq -c  | awk -v prop="\${threshold_pass_group}" '{if (\$1>=prop){print \$2}}' | cut -d "_" -f1,2 > positions_good_coverage.txt
    #Subset variant that have a support at least of X. If one sample has this support we keep the variant
    echo -e "Group-level filtering: Variant support ... ";
    awk -v OFS="\\t" -v FS="\\t" -v t_nv="${threshold_good_variant_support}" 'NR == FNR {  a[\$0]; next }{if(\$3"_"\$4 in a){sumval=\$9+\$10; if(sumval >= t_nv){print \$2}}}' positions_good_coverage.txt df_raw_variants_group_alt.tsv  | sort -u > variants_to_use_passed.txt
    #Subset the chosen variants 
    awk -v OFS="\\t" -v FS="\\t" 'NR == FNR {  a[\$0]; next }{if(\$2 in a){ print \$0}}' variants_to_use_passed.txt df_raw_variants.tsv | awk -v OFS="\\t" '{sumval=\$9+\$10;print \$1,\$2,sumval}' > df_nv.tsv
    
    #Return a table with depth too per sample
    cat df_nv.tsv | cut -f2 | cut -d "_" -f1,2 | sort -u > pos_nv.txt
    cut -f1,3,4,20 df_raw_variants_group.tsv | awk -v OFS="\\t" '{print \$1,\$2"_"\$3,\$4}'| sort -u > df_depths.tsv
    
    awk -v OFS="\\t" -v FS="\\t" 'NR == FNR {  a[\$0]; next }{if(\$2 in a){ print \$0}}' pos_nv.txt df_depths.tsv > df_nr.tsv
    ##############################
    ##### Matrix creation ########
    ##############################
    echo -e "Creating matrices ... ";
    Rscript /usr/local/bin/rscript_2.create_tabnr_tabnv.R
    echo -e "Subsetting matrices ... ";
    head -n1 Tab_NR.tsv > Mat_NR_${group}_${chr}.tsv
    tail -n +2 Tab_NR.tsv > body.tsv
    cat ${chosen_variants} | grep "^${chr}_" > mvariants.txt
    awk -v OFS="\\t" -v FS="\\t" 'NR == FNR {  a[\$0]; next }{if(\$1 in a){ print \$0}}' mvariants.txt body.tsv >> Mat_NR_${group}_${chr}.tsv
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