nextflow.enable.dsl=2

// IMPORT MODULES

include { PREPROCESS_VCF } from '../../modules/bcftools/filter_norm_addqual/main.nf' addParams( timestamp: params.timestamp )
include { MERGE_PROCESSED_VCF } from '../../modules/bcftools/merge_processed_vcf/main.nf' addParams( timestamp: params.timestamp )
include { CUSTOM_BAM_GROUP_PILEUP } from '../../modules/bioskryb/custom_bam_group_pileup/main.nf' addParams( timestamp: params.timestamp )
include { CREATE_TAB_NVNR } from '../../modules/bioskryb/create_tab_nvnr/main.nf' addParams( timestamp: params.timestamp )
include { SEQUOIA_BINOM_BETABINOM_TAB_NV_NR } from '../../modules/bioskryb/sequoia_binom_betabinom_tab_nv_nr/main.nf' addParams( timestamp: params.timestamp )
include { CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR } from '../../modules/bioskryb/concat_filter_binom_betabinom_tab_nv_nr/main.nf' addParams( timestamp: params.timestamp )
include { CUSTOM_RSCRIPT_SOMATICSNP_FILTER_1_SAMPLELEVEL_PROCESS_PILEUP_SAMPLE_CIGAR } from '../../modules/bioskryb/custom_rscript_somaticsnp_filter_1_samplelevel_process_pileup_sample_cigar/main.nf' addParams( timestamp: params.timestamp )
include { CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES } from '../../modules/bioskryb/custom_somatic_snpindel_filterrawtables/main.nf' addParams( timestamp: params.timestamp )
include { CUSTOM_CREATE_GROUP_LEVEL_TAB_DFS } from '../../modules/bioskryb/custom_create_group_level_tab_dfs/main.nf' addParams( timestamp: params.timestamp )
include { BULK_GET_VARIANTS_TO_FILTER } from '../../modules/bioskryb/custom_bulk_get_variants_to_filter/main.nf' addParams( timestamp: params.timestamp )
include { SEQUOIA } from '../../modules/sequoia/main.nf' addParams( timestamp: params.timestamp )
include { SUBSET_VCF_VARIANTS } from '../../modules/bioskryb/subset_vcf_variants/main.nf' addParams( timestamp: params.timestamp )
include { POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE } from '../../modules/bioskryb/custom_postprocess_sequoia_drawvafheat_tree/main.nf' addParams( timestamp: params.timestamp )
include { CUSTOM_VARIANT_FILTER_PROVENANCE } from '../../modules/bioskryb/custom_variant_filter_provenance/main.nf' addParams( timestamp: params.timestamp )
include { SUBSET_MERGED_VCF_CHOSEN_VARIANTS } from '../../modules/bioskryb/vep_chosen_variants/main.nf' addParams( timestamp: params.timestamp )
include { SPLIT_SUBSET_VCF_BY_CHR } from '../../modules/bioskryb/vep_chosen_variants/main.nf' addParams( timestamp: params.timestamp )
include { VEP_ANNOTATE } from '../../modules/bioskryb/vep_chosen_variants/main.nf' addParams( timestamp: params.timestamp )
include { SORT_INDEX_VEP_VCF } from '../../modules/bioskryb/vep_chosen_variants/main.nf' addParams( timestamp: params.timestamp )
include { MERGE_VEP_VCF_BY_GROUP } from '../../modules/bioskryb/vep_chosen_variants/main.nf' addParams( timestamp: params.timestamp )
include { FILTER_VEP_GERMLINE } from '../../modules/bioskryb/filter_vep_germline/main.nf' addParams( timestamp: params.timestamp )
include { CREATE_EMPTY_BULK_VARIANTS } from '../../modules/bioskryb/filter_chosen_variants_by_bulk/main.nf' addParams( timestamp: params.timestamp )
include { FILTER_CHOSEN_VARIANTS_BY_BULK } from '../../modules/bioskryb/filter_chosen_variants_by_bulk/main.nf' addParams( timestamp: params.timestamp )
include { GET_VARIANTS_FROM_MERGED_VCF } from '../../modules/bioskryb/get_variants_and_list_pos/main.nf' addParams( timestamp: params.timestamp )
include { GET_LIST_POS_FROM_CHOSEN_VARIANTS } from '../../modules/bioskryb/get_variants_and_list_pos/main.nf' addParams( timestamp: params.timestamp )
include { FILTER_DF_NV_BY_CHOSEN_VARIANTS } from '../../modules/bioskryb/get_variants_and_list_pos/main.nf' addParams( timestamp: params.timestamp )


workflow SOMATIC_SNP_INDEL_FILTERING_WF {

    Channel.fromPath( params.input_csv  ).
        splitCsv( header:true )
        .branch { row ->
            ch_bam: row.file_type == "bam"
                    return tuple( row.biosampleName, tuple ( row.file, row.file_index ), row.group  )
            ch_vcf: row.file_type == "vcf"
                    return tuple( row.biosampleName, tuple ( row.file ), row.group  )

        }
        .set { inputs }


    PREPROCESS_VCF (

        inputs.ch_vcf,
        params.reference,
        params.model_vcf,
        params.publish_dir,
        params.disable_publish

    )

    ch_input_merge_processed_vcf = PREPROCESS_VCF.out.vcf
    .map{ it -> [it[2],it[1]]}
    .groupTuple(by:0)
    .map{ it -> [it[0],it[1].flatten().collect()]}


    MERGE_PROCESSED_VCF (

        ch_input_merge_processed_vcf,
        params.reference,
        params.publish_dir,
        params.disable_publish

    )

    // Bulk filter right after merge: variants in bulk are removed; remainder drive pileup and downstream.
    if (params.bulk_vcf != "") {
        ch_bulk_vcf = Channel.fromPath(params.bulk_vcf)
        BULK_GET_VARIANTS_TO_FILTER (
            ch_bulk_vcf,
            params.reference,
            params.model_vcf,
            params.publish_dir,
            params.disable_publish
        )
        ch_bulk_variants_to_remove = BULK_GET_VARIANTS_TO_FILTER.out
    } else {
        CREATE_EMPTY_BULK_VARIANTS (params.publish_dir, params.enable_publish)
        ch_bulk_variants_to_remove = CREATE_EMPTY_BULK_VARIANTS.out.bulk_variants
    }

    // (group, path(merged_vcf)) -> (group, path(all_variants_${group}.txt))
    GET_VARIANTS_FROM_MERGED_VCF (
        MERGE_PROCESSED_VCF.out.merged_vcf,
        params.publish_dir,
        params.enable_publish
    )

    // (group, path(all_variants)), (path(bulk)) -> (group, path(all_variants), path(bulk))
    ch_filter_bulk_input_early = GET_VARIANTS_FROM_MERGED_VCF.out.all_variants
        .combine(ch_bulk_variants_to_remove)

    FILTER_CHOSEN_VARIANTS_BY_BULK (
        ch_filter_bulk_input_early,
        params.publish_dir,
        params.enable_publish
    )

    // (group, path(chosen_variants_filtered_bulk_*.txt))
    ch_chosen_after_bulk_early = FILTER_CHOSEN_VARIANTS_BY_BULK.out.chosen_variants

    // (group, path(chosen)) -> (group, path(list_pos_variant_*.txt))
    GET_LIST_POS_FROM_CHOSEN_VARIANTS (
        ch_chosen_after_bulk_early,
        params.publish_dir,
        params.enable_publish
    )

    ch_group_tables = PREPROCESS_VCF.out.query_table
    .map{ it -> [it[2],it[1]]}
    .groupTuple(by: 0)
    .map{ it -> [it[0],it[1].flatten().collect()]}

    ch_chr = Channel.of( params.chrs ).flatMap()

    // (group, sample_name, path(bam)), (group, path(list_pos)), (chr) -> (group, sample_name, path(bam), path(list_pos), chr)
    ch_input_bam_group_pileup = inputs.ch_bam
    .map{
        it -> [it[2],it[0],it[1]]
    }
    .combine(GET_LIST_POS_FROM_CHOSEN_VARIANTS.out.list_pos, by: 0)
    .combine(ch_chr)



    CUSTOM_BAM_GROUP_PILEUP (

        ch_input_bam_group_pileup,
        params.reference,
        params.publish_dir,
        params.disable_publish

    )

    ch_input_df_nr = CUSTOM_BAM_GROUP_PILEUP.out.df_nr
    .map{it ->
        [it[0],it[1],it[3]]
    }
    .groupTuple(by: [0,1])
    .map{it ->
        [it[0],it[1],it[2].flatten().collect()]
    }

    // (group, path(df_nv)), (group, path(chosen_after_bulk)) -> (group, path(df_nv), path(chosen))
    ch_input_filter_df_nv = MERGE_PROCESSED_VCF.out.df_nv
        .combine(ch_chosen_after_bulk_early, by: 0)

    FILTER_DF_NV_BY_CHOSEN_VARIANTS (
        ch_input_filter_df_nv,
        params.publish_dir,
        params.enable_publish
    )

    // (group, path(df_nv_filtered)), (group, chr, path(df_nr_files)) -> (group, chr, path(df_nv), path(df_nr_files))
    ch_input_create_tab_nvnr = FILTER_DF_NV_BY_CHOSEN_VARIANTS.out.df_nv_filtered
    .combine(ch_input_df_nr, by: 0)
    .map{ it ->
        [it[0], it[2], it[1], it[3].flatten().collect()]
    }


    CREATE_TAB_NVNR (

        ch_input_create_tab_nvnr,
        params.publish_dir,
        params.disable_publish

    )


    SEQUOIA_BINOM_BETABINOM_TAB_NV_NR(

        CREATE_TAB_NVNR.out,
        params.aggregated_min_mean_depth,
        params.aggregated_max_mean_depth,
        params.gender,
        params.publish_dir,
        params.disable_publish

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
        params.first_pass_binomial_cutoff,
        params.first_pass_betabinomial_cutoff,
        params.publish_dir,
        params.disable_publish

    )

    ch_input_vep = MERGE_PROCESSED_VCF.out.merged_vcf
        .combine(CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR.out.chosen_variants, by: 0)

    SUBSET_MERGED_VCF_CHOSEN_VARIANTS (
        ch_input_vep,
        params.publish_dir,
        params.enable_publish
    )

    ch_subset_by_chr = SUBSET_MERGED_VCF_CHOSEN_VARIANTS.out.subset_vcf
        .combine(Channel.of(params.chrs).flatMap())

    SPLIT_SUBSET_VCF_BY_CHR (
        ch_subset_by_chr,
        params.publish_dir,
        params.enable_publish
    )

    if (params.vep_cache_dir) {
        // Pass cache as path so Nextflow stages it from S3 into the container (required for Docker on EC2).
        // Must be the cache root (directory containing the species folder, e.g. .../VEP/).
        ch_vep_cache = file(params.vep_cache_dir, type: 'dir')
        VEP_ANNOTATE (
            SPLIT_SUBSET_VCF_BY_CHR.out.subset_vcf_chr,
            params.reference,
            params.vep_species,
            params.vep_assembly,
            ch_vep_cache,
            params.publish_dir,
            params.enable_publish
        )

        SORT_INDEX_VEP_VCF (
            VEP_ANNOTATE.out.vep_vcf_chr,
            params.publish_dir,
            params.enable_publish
        )

        // groupTuple(by: 0) on (group, chr, vcf, tbi) gives (group, [chrs], [vcfs], [tbis]); sort vcfs by chr order
        ch_vep_by_group = SORT_INDEX_VEP_VCF.out.vep_vcf_chr
            .groupTuple(by: 0)
            .map { group, list_chr, list_vcf, list_tbi -> [group, [list_chr, list_vcf].transpose().sort{ a, b -> a[0] <=> b[0] }.collect{ it[1] }] }

        MERGE_VEP_VCF_BY_GROUP (
            ch_vep_by_group,
            params.publish_dir,
            params.enable_publish
        )

        FILTER_VEP_GERMLINE (
            MERGE_VEP_VCF_BY_GROUP.out.vep_vcf,
            params.vep_max_af_1kg,
            params.vep_filter_by_existing_variation,
            params.publish_dir,
            params.enable_publish
        )
        ch_chosen_variants = FILTER_VEP_GERMLINE.out.chosen_variants
    } else {
        ch_chosen_variants = CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR.out.chosen_variants
    }

    // Bulk filter already applied early (after MERGE); ch_chosen_variants is CONCAT or VEP output.
    ch_input_rscript_filter_1 = CUSTOM_BAM_GROUP_PILEUP.out.pileup
    .combine(ch_chosen_variants, by: 0)

    CUSTOM_RSCRIPT_SOMATICSNP_FILTER_1_SAMPLELEVEL_PROCESS_PILEUP_SAMPLE_CIGAR (

        ch_input_rscript_filter_1,
        params.cutoff_mq_hq,
        params.cutoff_bq_hq,
        params.cutoff_bps_start,
        params.num_lines_read_pileup,
        params.read_length,
        params.publish_dir,
        params.enable_publish

    )

    ch_input_filter_tables = CUSTOM_RSCRIPT_SOMATICSNP_FILTER_1_SAMPLELEVEL_PROCESS_PILEUP_SAMPLE_CIGAR.out
    .groupTuple(by :[0,1])
    .map{it ->
        [it[0],it[1],it[2].flatten().collect()]
    }
    .combine(ch_chosen_variants, by: 0)

    CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES (

        ch_input_filter_tables,
        params.cutoff_as,
        params.cutoff_prop_clipped_reads,
        params.cutoff_prop_bp_under,
        params.cutoff_prop_bp_upper,
        params.cutoff_sd_indiv,
        params.cutoff_mad_indiv,
        params.cutoff_sd_both,
        params.cutoff_mad_both,
        params.cutoff_sd_extreme,
        params.cutoff_mad_extreme,
        params.cutoff_prop_cells_goodcov_group,
        params.cutoff_goodcov_depth,
        params.cutoff_numreads_variant_manual,
        params.cutoff_num_hq_fragments_forward,
        params.cutoff_num_hq_fragments_reverse,
        params.disable_qc,
        params.publish_dir,
        params.disable_publish

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
        params.publish_dir,
        params.enable_publish

    )

    SEQUOIA (

        CUSTOM_CREATE_GROUP_LEVEL_TAB_DFS.out.tabs,
        params.reference,
        params.second_pass_binomial_cutoff,
        params.second_pass_betabinomial_cutoff_rho_snp,
        params.second_pass_betabinomial_cutoff_rho_indel,
        params.aggregated_hq_min_mean_depth,
        params.aggregated_hq_max_mean_depth,
        params.publish_dir,
        params.enable_publish

    )


    ch_input_vcf_subset = PREPROCESS_VCF.out.vcf
    .map{it ->

        [it[2],it[0],it[1]]

    }.combine(SEQUOIA.out.df,by:0)

    SUBSET_VCF_VARIANTS (

        ch_input_vcf_subset,
        params.reference,
        params.publish_dir,
        params.enable_publish

    )

    ch_all_df_gt = PREPROCESS_VCF.out.df_gt
    .map{
        it -> [it[1]]
    }
    .collect()

    POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE (

        SEQUOIA.out.bundle_post_vaf_tree,
        ch_all_df_gt,
        params.publish_dir,
        params.enable_publish

    )


    ch_df_gt = PREPROCESS_VCF.out.df_gt
    .map{
        it -> [it[2],it[1]]
    }
    .groupTuple(by:0)
    .map{
        it -> [it[0],it[1].flatten().collect()]
    }


    ch_qc_tables = CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES.out.df_pass
    .map{
        it -> [it[0],it[2]]
    }
    .groupTuple(by:0)
    .map{
        it -> [it[0],it[1].flatten().collect()]
    }

    // Provenance: (group, res_df, chosen_after_binomial, chosen_before_bulk, chosen_after_bulk, chosen_after_vep, path(Mat_NV))
    ch_input_variant_provenance = CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR.out.res_df
    .combine(CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR.out.chosen_variants, by: 0)
    .combine(GET_VARIANTS_FROM_MERGED_VCF.out.all_variants, by: 0)
    .combine(ch_chosen_after_bulk_early, by: 0)
    .combine(ch_chosen_variants, by: 0)
    .combine(CUSTOM_CREATE_GROUP_LEVEL_TAB_DFS.out.tabs.map { it -> [it[0], it[1]] }, by: 0)

    CUSTOM_VARIANT_FILTER_PROVENANCE (

        ch_input_variant_provenance,
        params.publish_dir,
        params.enable_publish

    )

    ch_all_subset_vcf = SUBSET_VCF_VARIANTS.out.collect()


}
