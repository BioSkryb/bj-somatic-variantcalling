nextflow.enable.dsl=2

// IMPORT MODULES
include { PREPROCESS_VCF } from '../../modules/bcftools/filter_norm_addqual/main.nf'
include { MERGE_PROCESSED_VCF } from '../../modules/bcftools/merge_processed_vcf/main.nf'
include { CUSTOM_GET_LIST_POS_GROUP } from '../../modules/bioskryb/custom_get_list_pos_group/main.nf'
include { CUSTOM_BAM_GROUP_PILEUP } from '../../modules/bioskryb/custom_bam_group_pileup/main.nf'
include { CREATE_TAB_NVNR } from '../../modules/bioskryb/create_tab_nvnr/main.nf'
include { SEQUOIA_BINOM_BETABINOM_TAB_NV_NR } from '../../modules/bioskryb/sequoia_binom_betabinom_tab_nv_nr/main.nf'
include { CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR } from '../../modules/bioskryb/concat_filter_binom_betabinom_tab_nv_nr/main.nf'
include { CUSTOM_CREATE_GROUP_LEVEL_TAB_DFS } from '../../modules/bioskryb/custom_create_group_level_tab_dfs/main.nf'
include { CUSTOM_SPLIT_BAM_CHR } from '../../modules/bioskryb/custom_split_bam_chr/main.nf'
include { CUSTOM_BAM_GROUP_PILEUP_CHR } from '../../modules/bioskryb/custom_bam_group_pileup_chr/main.nf'
include { CUSTOM_RSCRIPT_SOMATICSNP_FILTER_1_SAMPLELEVEL_PROCESS_PILEUP_SAMPLE_CIGAR } from '../../modules/bioskryb/custom_rscript_somaticsnp_filter_1_samplelevel_process_pileup_sample_cigar/main.nf'
include { CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES } from '../../modules/bioskryb/custom_somatic_snpindel_filterrawtables/main.nf'
include { SEQUOIA } from '../../modules/sequoia/main.nf'
include { SUBSET_VCF_VARIANTS } from '../../modules/bioskryb/subset_vcf_variants/main.nf'
include { POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE } from '../../modules/bioskryb/custom_postprocess_sequoia_drawvafheat_tree/main.nf'

workflow SOMATIC_VARIANT_WORKFLOW_Heuristic_Filter {
    take:
    ch_bam
    ch_vcf
    reference
    chrs
    cutoff_mq_hq
    cutoff_bq_hq
    cutoff_bps_start
    cutoff_as
    cutoff_prop_clipped_reads
    cutoff_prop_bp_under
    cutoff_sd_indiv
    cutoff_mad_indiv
    cutoff_sd_both
    cutoff_mad_both
    cutoff_sd_extreme
    cutoff_mad_extreme
    cutoff_prop_cells_goodcov_group
    cutoff_goodcov_depth
    cutoff_numreads_variant_manual
    ch_model_vcf
    first_pass_binomial_cutoff
    first_pass_betabinomial_cutoff
    num_lines_read_pileup
    second_pass_binomial_cutoff
    second_pass_betabinomial_cutoff_rho_snp
    second_pass_betabinomial_cutoff_rho_indel
    aggregated_hq_min_mean_depth
    aggregated_hq_max_mean_depth
    aggregated_min_mean_depth
    aggregated_max_mean_depth
    gender
    publish_dir
    disable_publish
    enable_publish

    main:

    PREPROCESS_VCF(
        ch_vcf,
        reference,
        publish_dir,
        ch_model_vcf,
        disable_publish
    )

    ch_input_merge_processed_vcf = PREPROCESS_VCF.out.vcf
    .map{ it -> [it[2],it[1]]}
    .groupTuple(by:0)
    .map{ it -> [it[0],it[1].flatten().collect()]}


    MERGE_PROCESSED_VCF (

        ch_input_merge_processed_vcf,
        reference,
        publish_dir,
        disable_publish 

    )

    ch_group_tables = PREPROCESS_VCF.out.query_table
    .map{ it -> [it[2],it[1]]}
    .groupTuple(by: 0)
    .map{ it -> [it[0],it[1].flatten().collect()]}

    // ch_group_tables.view()

    CUSTOM_GET_LIST_POS_GROUP (
        ch_group_tables,
        publish_dir,
        disable_publish
    )

    ch_chr = Channel.of( chrs ).flatMap()

    ch_input_bam_group_pileup = ch_bam
    .map{
        it -> [it[2],it[0],it[1]]
    }
    .combine(CUSTOM_GET_LIST_POS_GROUP.out,by:0)
    .combine(ch_chr)

    CUSTOM_BAM_GROUP_PILEUP (

        ch_input_bam_group_pileup,
        reference,
        publish_dir,
        disable_publish

    )

    ch_input_df_nr = CUSTOM_BAM_GROUP_PILEUP.out.df_nr
    .map{it ->
        [it[0],it[1],it[3]]
    }
    .groupTuple(by: [0,1])
    .map{it ->
        [it[0],it[1],it[2].flatten().collect()]
    }

    ch_input_create_tab_nvnr = MERGE_PROCESSED_VCF.out.df_nv
    .combine(ch_input_df_nr,by:0)
    .map{it ->
        [it[0],it[2],it[1],it[3].flatten().collect()]
    }

    CREATE_TAB_NVNR (

        ch_input_create_tab_nvnr,
        publish_dir,
        disable_publish

    )

    SEQUOIA_BINOM_BETABINOM_TAB_NV_NR(

        CREATE_TAB_NVNR.out,
        aggregated_min_mean_depth,
        aggregated_max_mean_depth,
        gender,
        publish_dir,
        disable_publish

    )

    ch_input_concat_filter_bb = SEQUOIA_BINOM_BETABINOM_TAB_NV_NR.out.df_filter
    .map{ it->
        [it[0],it[2]]
    }
    .groupTuple(by:[0])
    .map{ it->
        [it[0],it[1].flatten().collect()]
    }

    CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR (

        ch_input_concat_filter_bb,
        first_pass_binomial_cutoff,
        first_pass_betabinomial_cutoff,
        publish_dir,
        disable_publish

    )

    ch_input_rscript_filter_1 = CUSTOM_BAM_GROUP_PILEUP.out.pileup
    .combine(CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR.out.chosen_variants,by: 0)



    CUSTOM_RSCRIPT_SOMATICSNP_FILTER_1_SAMPLELEVEL_PROCESS_PILEUP_SAMPLE_CIGAR (

        ch_input_rscript_filter_1,
        cutoff_mq_hq,
        cutoff_bq_hq,
        cutoff_bps_start,
        num_lines_read_pileup,
        publish_dir,
        disable_publish

    )

    ch_input_filter_tables = CUSTOM_RSCRIPT_SOMATICSNP_FILTER_1_SAMPLELEVEL_PROCESS_PILEUP_SAMPLE_CIGAR.out
    .groupTuple(by :[0,1])
    .map{it ->
        [it[0],it[1],it[2].flatten().collect()]
    }
    .combine(CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR.out.chosen_variants,by: 0)

    CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES (
        ch_input_filter_tables,
        cutoff_as,
        cutoff_prop_clipped_reads,
        cutoff_prop_bp_under,
        cutoff_sd_indiv,
        cutoff_mad_indiv,
        cutoff_sd_both,
        cutoff_mad_both,
        cutoff_sd_extreme,
        cutoff_mad_extreme,
        cutoff_prop_cells_goodcov_group,
        cutoff_goodcov_depth,
        cutoff_numreads_variant_manual,
        publish_dir,
        enable_publish
    )

    ch_input_tabs_group_level = CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES.out.tabs
    .map{ it->
        [it[0],it[2],it[3]]
    }
    .groupTuple(by:0)
    .map{ it ->
        [it[0],it[1].flatten().collect(),it[2].flatten().collect()]
    }

    CUSTOM_CREATE_GROUP_LEVEL_TAB_DFS (

        ch_input_tabs_group_level,
        publish_dir,
        enable_publish

    )

    SEQUOIA (
        CUSTOM_CREATE_GROUP_LEVEL_TAB_DFS.out.tabs,
        reference,
        second_pass_binomial_cutoff,
        second_pass_betabinomial_cutoff_rho_snp,
        second_pass_betabinomial_cutoff_rho_indel,
        aggregated_hq_min_mean_depth,
        aggregated_hq_max_mean_depth,
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

    POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE (

        SEQUOIA.out.bundle_post_vaf_tree,
        publish_dir,
        enable_publish

    )

    emit:
    vcf = SUBSET_VCF_VARIANTS.out.vcf
    sequoia_results = SEQUOIA.out.df
}