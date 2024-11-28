nextflow.enable.dsl=2

// IMPORT MODULES
include { PREPROCESS_VCF } from '../../modules/bcftools/filter_norm_addqual/main.nf' addParams( timestamp: params.timestamp )
include { CUSTOM_SPLIT_QUERY_TABLE_CHR } from '../../modules/bioskryb/custom_split_query_table_chr/main.nf' addParams( timestamp: params.timestamp )
include { CUSTOM_SPLIT_BAM_CHR } from '../../modules/bioskryb/custom_split_bam_chr/main.nf' addParams( timestamp: params.timestamp )
include { CUSTOM_BAM_PILEUP_FILTER } from '../../modules/bioskryb/custom_bam_pileup_filter/main.nf' addParams( timestamp: params.timestamp )
include { CUSTOM_RSCRIPT_SOMATICSNP_FILTER_1_CELL_LEVEL_CREATE_TABLE } from '../../modules/bioskryb/custom_rscript_somaticsnp_filter_1_celllevel_create_table/main.nf' addParams( timestamp: params.timestamp )
include { CUSTOM_RSCRIPT_SOMATICSNP_FILTER_2_CELL_LEVEL_CONCAT_FILTER_TABLE } from '../../modules/bioskryb/custom_rscript_somaticsnp_filter_2_celllevel_concat_filter_table/main.nf' addParams( timestamp: params.timestamp )
include { CUSTOM_GET_LIST_POS_GROUP } from '../../modules/bioskryb/custom_get_list_pos_group/main.nf' addParams( timestamp: params.timestamp )
include { CUSTOM_BAM_GROUP_PILEUP_CHR } from '../../modules/bioskryb/custom_bam_group_pileup_chr/main.nf' addParams( timestamp: params.timestamp )
include { CUSTOM_RSCRIPT_SOMATICSNP_FILTER_3_GROUPLEVEL_PROCESS_PILEUP_SAMPLE } from '../../modules/bioskryb/custom_rscript_somaticsnp_filter_3_grouplevel_process_pileup_sample/main.nf' addParams( timestamp: params.timestamp )
include { CUSTOM_RSCRIPT_SOMATICSNP_FILTER_4_CREATE_TABNR_TABNV } from '../../modules/bioskryb/custom_rscript_somaticsnp_filter_4_create_tabnr_tabnv/main.nf' addParams( timestamp: params.timestamp )
include { SEQUOIA } from '../../modules/sequoia/main.nf' addParams( timestamp: params.timestamp )
include { SUBSET_VCF_VARIANTS } from '../../modules/bioskryb/subset_vcf_variants/main.nf' addParams( timestamp: params.timestamp )

workflow SOMATIC_VARIANT_WORKFLOW_Heuristic_Filter {
    take:
    ch_bam
    ch_vcf
    reference
    chrs
    num_lines_split_query_table
    cutoff_as
    cutoff_prop_clipped_reads
    cutoff_num_hq_frag
    filter_bp_pos
    cutoff_prop_bp
    cutoff_sd_indiv
    cutoff_mad_indiv
    cutoff_sd_both
    cutoff_mad_both
    cutoff_sd_extreme
    cutoff_mad_extreme
    num_lines_get_list_pos_group
    cutoff_depth_manual
    cutoff_numreads_variant_manual
    cutoff_prev_na_manual
    cutoff_prev_var_manual
    sequoia_cutoff_binomial
    sequoia_cutoff_rho
    publish_dir
    disable_publish
    enable_publish

    main:

    PREPROCESS_VCF(
        ch_vcf,
        reference,
        publish_dir,
        disable_publish
    )
    
    ch_query_table = PREPROCESS_VCF.out.query_table

    ch_chr = Channel.of( chrs ).flatMap()

    ch_query_table_chr = ch_query_table
        .combine( ch_chr )

    ch_bam_chr = ch_bam
        .combine( ch_chr )
    
    CUSTOM_SPLIT_QUERY_TABLE_CHR(
        ch_query_table_chr,
        num_lines_split_query_table,
        publish_dir,
        disable_publish
    )
    
    CUSTOM_SPLIT_BAM_CHR(
        ch_bam_chr,
        publish_dir,
        disable_publish
    )
    
    ch_query_table = CUSTOM_SPLIT_QUERY_TABLE_CHR.out
        .flatten()
        .map { it -> [ it.BaseName.replaceFirst(/.*_sample_/, "").replaceFirst(/_chr_.*/,""), it.BaseName.replaceFirst(/.*_chr_/, "").replaceFirst(/_file_.*/,"") , it.BaseName.replaceFirst(/.*_file_/, "").replaceFirst(/\.txt/,"") , it ] }

    ch_bam_chr_ipile = CUSTOM_SPLIT_BAM_CHR.out.bam

    ch_input_pileup = ch_query_table
        .combine(ch_bam_chr_ipile, by: [0,1] )
    
    CUSTOM_BAM_PILEUP_FILTER(
        ch_input_pileup,
        reference,
        publish_dir,
        disable_publish
    )
    
    ch_bam_pileup = CUSTOM_BAM_PILEUP_FILTER.out
    
    ch_input_rscript = ch_bam_pileup
        .combine(CUSTOM_SPLIT_BAM_CHR.out.df, by: [0,1] )

    CUSTOM_RSCRIPT_SOMATICSNP_FILTER_1_CELL_LEVEL_CREATE_TABLE(
        ch_input_rscript,
        publish_dir,
        disable_publish
    )
    
    ch_concat_tables = CUSTOM_RSCRIPT_SOMATICSNP_FILTER_1_CELL_LEVEL_CREATE_TABLE.out
        .groupTuple(by: [0,1])
    
    CUSTOM_RSCRIPT_SOMATICSNP_FILTER_2_CELL_LEVEL_CONCAT_FILTER_TABLE(
        ch_concat_tables,
        cutoff_as,
        cutoff_prop_clipped_reads,
        cutoff_num_hq_frag,
        filter_bp_pos,
        cutoff_prop_bp,
        cutoff_sd_indiv,
        cutoff_mad_indiv,
        cutoff_sd_both,
        cutoff_mad_both,
        cutoff_sd_extreme,
        cutoff_mad_extreme,
        publish_dir,
        enable_publish
    )
    
    ch_group_level_samples_chr = CUSTOM_RSCRIPT_SOMATICSNP_FILTER_2_CELL_LEVEL_CONCAT_FILTER_TABLE.out.res_filter
        .map{ it ->
            [it[0],it[2]]
        }
        .groupTuple(by: 0)
        .combine( ch_chr )

    CUSTOM_GET_LIST_POS_GROUP(
        ch_group_level_samples_chr,
        num_lines_get_list_pos_group,
        publish_dir,
        disable_publish
    )
    
    ch_bam_chr_shuf = ch_bam_chr
        .map { it ->
            [it[2],it[3],it[0],it[1]]
        }
    
    ch_list_out = CUSTOM_GET_LIST_POS_GROUP.out
        .flatten()
        .map { it -> [ it.BaseName.replaceFirst(/.*_group_/, "").replaceFirst(/_chr_.*/,""), it.BaseName.replaceFirst(/.*_chr_/, "").replaceFirst(/_file_.*/,"") , it.BaseName.replaceFirst(/.*_file_/, "").replaceFirst(/\.txt/,"") , it ] }
    
    ch_input_bam_list_group = ch_bam_chr_shuf
        .combine(ch_list_out,by: [0,1] )

    CUSTOM_BAM_GROUP_PILEUP_CHR(
        ch_input_bam_list_group,
        reference,
        publish_dir,
        disable_publish
    )

    CUSTOM_RSCRIPT_SOMATICSNP_FILTER_3_GROUPLEVEL_PROCESS_PILEUP_SAMPLE(
        CUSTOM_BAM_GROUP_PILEUP_CHR.out,
        publish_dir,
        enable_publish
    )
    
    ch_group_input_create_tabs_sequoia = CUSTOM_RSCRIPT_SOMATICSNP_FILTER_3_GROUPLEVEL_PROCESS_PILEUP_SAMPLE.out
        .groupTuple(by: [0,1])

    CUSTOM_RSCRIPT_SOMATICSNP_FILTER_4_CREATE_TABNR_TABNV(
        ch_group_input_create_tabs_sequoia,
        cutoff_depth_manual,
        cutoff_numreads_variant_manual,
        cutoff_prev_na_manual,
        cutoff_prev_var_manual,
        publish_dir,
        enable_publish
    )
    
    ch_input_sequoia = CUSTOM_RSCRIPT_SOMATICSNP_FILTER_4_CREATE_TABNR_TABNV.out.clean
        .groupTuple(by:0)
        .map{ it->
            [it[0],it[1].flatten().collect()]
        }
    
    SEQUOIA(
        ch_input_sequoia,
        reference,
        sequoia_cutoff_binomial,
        sequoia_cutoff_rho,
        publish_dir,
        enable_publish
    )
    
    ch_input_vcf_subset = PREPROCESS_VCF.out.vcf
        .map{it ->
            [it[2],it[0],it[1]]
        }.combine(SEQUOIA.out.df,by:0)
    
    SUBSET_VCF_VARIANTS(
        ch_input_vcf_subset,
        reference,
        publish_dir,
        enable_publish
    )

    emit:
    vcf = SUBSET_VCF_VARIANTS.out.vcf
    sequoia_results = SEQUOIA.out.df
}