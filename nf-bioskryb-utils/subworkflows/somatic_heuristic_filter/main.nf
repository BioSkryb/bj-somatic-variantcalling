nextflow.enable.dsl=2

// ============================================================================
// IMPORT MODULES
// ============================================================================

include { PREPROCESS_VCF                                                         } from '../../modules/bcftools/filter_norm_addqual/main.nf'                                                              addParams( timestamp: params.timestamp )
include { MERGE_PROCESSED_VCF                                                    } from '../../modules/bcftools/merge_processed_vcf/main.nf'                                                              addParams( timestamp: params.timestamp )
include { CUSTOM_BAM_GROUP_PILEUP                                                } from '../../modules/bioskryb/custom_bam_group_pileup/main.nf'                                                          addParams( timestamp: params.timestamp )
include { CREATE_TAB_NVNR                                                        } from '../../modules/bioskryb/create_tab_nvnr/main.nf'                                                                  addParams( timestamp: params.timestamp )
include { SEQUOIA_BINOM_BETABINOM_TAB_NV_NR                                     } from '../../modules/bioskryb/sequoia_binom_betabinom_tab_nv_nr/main.nf'                                                addParams( timestamp: params.timestamp )
include { CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR                               } from '../../modules/bioskryb/concat_filter_binom_betabinom_tab_nv_nr/main.nf'                                          addParams( timestamp: params.timestamp )
include { CUSTOM_RSCRIPT_SOMATICSNP_FILTER_1_SAMPLELEVEL_PROCESS_PILEUP_SAMPLE_CIGAR } from '../../modules/bioskryb/custom_rscript_somaticsnp_filter_1_samplelevel_process_pileup_sample_cigar/main.nf' addParams( timestamp: params.timestamp )
include { CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES                               } from '../../modules/bioskryb/custom_somatic_snpindel_filterrawtables/main.nf'                                          addParams( timestamp: params.timestamp )
include { CUSTOM_CREATE_GROUP_LEVEL_TAB_DFS                                     } from '../../modules/bioskryb/custom_create_group_level_tab_dfs/main.nf'                                                addParams( timestamp: params.timestamp )
include { BULK_GET_VARIANTS_TO_FILTER                                            } from '../../modules/bioskryb/custom_bulk_get_variants_to_filter/main.nf'                                               addParams( timestamp: params.timestamp )
include { SEQUOIA_SECOND_FILTER                                                  } from '../../modules/sequoia/main.nf'                                                                                   addParams( timestamp: params.timestamp )
include { SUBSET_VCF_VARIANTS                                                    } from '../../modules/bioskryb/subset_vcf_variants/main.nf'                                                              addParams( timestamp: params.timestamp )
include { POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE                                   } from '../../modules/bioskryb/custom_postprocess_sequoia_drawvafheat_tree/main.nf'                                      addParams( timestamp: params.timestamp )
include { CUSTOM_VARIANT_FILTER_PROVENANCE                                       } from '../../modules/bioskryb/custom_variant_filter_provenance/main.nf'                                                 addParams( timestamp: params.timestamp )
include { LIST_SAMPLES_FROM_GROUP_VCF                                            } from '../../modules/bioskryb/list_samples_from_vcf/main.nf'                                                            addParams( timestamp: params.timestamp )
include { ANNOTATE_SAMPLE_VCF                                                    } from '../../modules/bioskryb/annotate_sample_vcf/main.nf'                                                              addParams( timestamp: params.timestamp )
include { GENOTYPE_TABLE_FROM_ANNOTATED_VCF                                      } from '../../modules/bioskryb/genotype_table_from_annotated_vcf/main.nf'                                               addParams( timestamp: params.timestamp )
include { CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF                               } from '../../modules/bioskryb/create_nr_nv_matrices_from_annotated_vcf/main.nf'                                        addParams( timestamp: params.timestamp )
include { SEQUOIA_PHYLOGENY_SNV                                                  } from '../../modules/bioskryb/sequoia_phylogeny_snv/main.nf'                                                            addParams( timestamp: params.timestamp )
include { SEQUOIA_PHYLOGENY_INDEL                                                } from '../../modules/bioskryb/sequoia_phylogeny_indel/main.nf'                                                          addParams( timestamp: params.timestamp )
include { SEQUOIA_PHYLOGENY_BOTH                                                 } from '../../modules/bioskryb/sequoia_phylogeny_both/main.nf'                                                           addParams( timestamp: params.timestamp )
include { TREES_COMPARE_SIMILARITIES                                             } from '../../modules/bioskryb/trees_compare_similarities/main.nf'                                                       addParams( timestamp: params.timestamp )
include { SEQUOIA_VARIANT_PLACEMENT_SNV                                          } from '../../modules/bioskryb/sequoia_variant_placement_snv/main.nf'                                                    addParams( timestamp: params.timestamp )
include { SEQUOIA_VARIANT_PLACEMENT_INDEL                                        } from '../../modules/bioskryb/sequoia_variant_placement_indel/main.nf'                                                  addParams( timestamp: params.timestamp )
include { SEQUOIA_VARIANT_PLACEMENT_BOTH                                         } from '../../modules/bioskryb/sequoia_variant_placement_both/main.nf'                                                   addParams( timestamp: params.timestamp )
include { SUBSET_ANNOTATED_VCFS_FOR_MUTSIG                                       } from '../../modules/bioskryb/subset_annotated_vcfs_for_mutsig/main.nf'                                                addParams( timestamp: params.timestamp )
include { SIGPROFILER_ASSIGNMENT                                                  } from '../../modules/bioskryb/sigprofiler_assignment/main.nf'                                                          addParams( timestamp: params.timestamp )
include { MERGE_SIGNATURE_ACTIVITIES                                             } from '../../modules/bioskryb/merge_signature_activities/main.nf'                                                       addParams( timestamp: params.timestamp )
include { PLOT_ASSIGNED_SIGNATURE_ACTIVITIES as PLOT_ZERO_FILTERED_SIGNATURE_ACTIVITIES   } from '../../modules/bioskryb/plot_assigned_signature_activities/main.nf' addParams( timestamp: params.timestamp )
include { PLOT_ASSIGNED_SIGNATURE_ACTIVITIES as PLOT_COSINE_FILTERED_SIGNATURE_ACTIVITIES } from '../../modules/bioskryb/plot_assigned_signature_activities/main.nf' addParams( timestamp: params.timestamp )
include { FILTER_VARIANTS_BY_SIG_PROBABILITY                                     } from '../../modules/bioskryb/filter_variants_by_sig_probability/main.nf'                                              addParams( timestamp: params.timestamp )
include { MERGE_SIG_COSINE_SIMILARITIES                                          } from '../../modules/bioskryb/merge_sig_cosine_similarities/main.nf'                                                    addParams( timestamp: params.timestamp )
include { FILTER_ACTIVITIES_BY_COSINE                                            } from '../../modules/bioskryb/filter_activities_by_cosine/main.nf'                                                      addParams( timestamp: params.timestamp )
include { SUBSET_MERGED_VCF_CHOSEN_VARIANTS                                      } from '../../modules/bioskryb/vep_chosen_variants/main.nf'                                                             addParams( timestamp: params.timestamp )
include { SPLIT_SUBSET_VCF_BY_CHR                                                } from '../../modules/bioskryb/vep_chosen_variants/main.nf'                                                             addParams( timestamp: params.timestamp )
include { VEP_ANNOTATE                                                           } from '../../modules/bioskryb/vep_chosen_variants/main.nf'                                                              addParams( timestamp: params.timestamp )
include { SORT_INDEX_VEP_VCF                                                     } from '../../modules/bioskryb/vep_chosen_variants/main.nf'                                                             addParams( timestamp: params.timestamp )
include { MERGE_VEP_VCF_BY_GROUP                                                 } from '../../modules/bioskryb/vep_chosen_variants/main.nf'                                                             addParams( timestamp: params.timestamp )
include { FILTER_VEP_GERMLINE                                                    } from '../../modules/bioskryb/filter_vep_germline/main.nf'                                                              addParams( timestamp: params.timestamp )
include { CREATE_EMPTY_BULK_VARIANTS                                             } from '../../modules/bioskryb/filter_chosen_variants_by_bulk/main.nf'                                                   addParams( timestamp: params.timestamp )
include { FILTER_CHOSEN_VARIANTS_BY_BULK                                         } from '../../modules/bioskryb/filter_chosen_variants_by_bulk/main.nf'                                                   addParams( timestamp: params.timestamp )
include { GET_VARIANTS_FROM_MERGED_VCF                                           } from '../../modules/bioskryb/get_variants_and_list_pos/main.nf'                                                        addParams( timestamp: params.timestamp )
include { GET_LIST_POS_FROM_CHOSEN_VARIANTS                                      } from '../../modules/bioskryb/get_variants_and_list_pos/main.nf'                                                        addParams( timestamp: params.timestamp )
include { FILTER_DF_NV_BY_CHOSEN_VARIANTS                                        } from '../../modules/bioskryb/get_variants_and_list_pos/main.nf'                                                        addParams( timestamp: params.timestamp )
include { COMPILE_MASTER_REPORT                                                  } from '../../modules/bioskryb/compile_master_report/main.nf'                                                            addParams( timestamp: params.timestamp )
include { IDENTIFY_GERMLINE_FROM_STATS                                           } from '../../modules/bioskryb/identify_germline_from_stats/main.nf'                                                          addParams( timestamp: params.timestamp )
include { EXTRACT_GERMLINE_PREVALENCE_TABLE                                      } from '../../modules/bioskryb/extract_germline_prevalence_table/main.nf'                                                     addParams( timestamp: params.timestamp )
include { PLOT_GERMLINE_PREVALENCE_DISTRIBUTIONS                                 } from '../../modules/bioskryb/plot_germline_prevalence_distributions/main.nf'                                               addParams( timestamp: params.timestamp )
include { SUBSET_MERGED_VCF_HIGH_CONFIDENCE_GERMLINE_FROM_STATS                  } from '../../modules/bioskryb/subset_merged_vcf_high_confidence_germline_from_stats/main.nf'                               addParams( timestamp: params.timestamp )
include { CREATE_ADO_TABLE_FROM_GERMLINE_VCF                                     } from '../../modules/bioskryb/create_ado_table_from_germline_vcf/main.nf'                                                   addParams( timestamp: params.timestamp )
include { SUMMARIZE_ADO_INTERVALS                                                 } from '../../modules/bioskryb/ado/summarize_ado_intervals_r/main.nf'                                                        addParams( timestamp: params.timestamp )
include { CONCAT_SUMMARY_ADO_INTERVALS_LABELED as CONCAT_ADO_STATS               } from '../../modules/bioskryb/concat_summary_ado_intervals_labeled/main.nf'                                                addParams( timestamp: params.timestamp )
include { CONCAT_SUMMARY_ADO_INTERVALS_LABELED as CONCAT_ADO_VEP                 } from '../../modules/bioskryb/concat_summary_ado_intervals_labeled/main.nf'                                                addParams( timestamp: params.timestamp )
include { CONCAT_SUMMARY_ADO_INTERVALS_LABELED as CONCAT_ADO_BULK                } from '../../modules/bioskryb/concat_summary_ado_intervals_labeled/main.nf'                                                addParams( timestamp: params.timestamp )
include { PLOT_ADO_GERMLINE_COMPARISON                                           } from '../../modules/bioskryb/plot_ado_germline_comparison/main.nf'                                                          addParams( timestamp: params.timestamp )

// ============================================================================
// WORKFLOW
// ============================================================================

workflow SOMATIC_SNP_INDEL_FILTERING_WF {

    Channel.fromPath( params.input_csv ).
        splitCsv( header:true )
        .multiMap { row ->
            ch_bam: tuple( row.biosampleName, tuple( row.bam, row.bam + ".bai" ), row.groups )
            ch_vcf: tuple( row.biosampleName, tuple( row.vcf ), row.groups )
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
        .map{ it -> [it[2],it[1]] }
        .groupTuple(by:0)
        .map{ it -> [it[0],it[1].flatten().collect()] }

    MERGE_PROCESSED_VCF (
        ch_input_merge_processed_vcf,
        params.reference,
        params.publish_dir,
        params.disable_publish
    )

    // Bulk filter right after merge
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

    GET_VARIANTS_FROM_MERGED_VCF (
        MERGE_PROCESSED_VCF.out.merged_vcf,
        params.publish_dir,
        params.enable_publish
    )

    ch_filter_bulk_input_early = GET_VARIANTS_FROM_MERGED_VCF.out.all_variants
        .combine(ch_bulk_variants_to_remove)

    FILTER_CHOSEN_VARIANTS_BY_BULK (
        ch_filter_bulk_input_early,
        params.publish_dir,
        params.enable_publish
    )

    ch_chosen_after_bulk_early = FILTER_CHOSEN_VARIANTS_BY_BULK.out.chosen_variants

    GET_LIST_POS_FROM_CHOSEN_VARIANTS (
        ch_chosen_after_bulk_early,
        params.publish_dir,
        params.enable_publish
    )

    ch_group_tables = PREPROCESS_VCF.out.query_table
        .map{ it -> [it[2],it[1]] }
        .groupTuple(by: 0)
        .map{ it -> [it[0],it[1].flatten().collect()] }

    ch_chr = Channel.of( params.chrs instanceof List ? params.chrs : params.chrs.tokenize(',') ).flatMap()

    ch_input_bam_group_pileup = inputs.ch_bam
        .map{ it -> [it[2],it[0],it[1]] }
        .combine(GET_LIST_POS_FROM_CHOSEN_VARIANTS.out.list_pos, by: 0)
        .combine(ch_chr)

    CUSTOM_BAM_GROUP_PILEUP (
        ch_input_bam_group_pileup,
        params.reference,
        params.publish_dir,
        params.disable_publish
    )

    ch_input_df_nr = CUSTOM_BAM_GROUP_PILEUP.out.df_nr
        .map{ it -> [it[0],it[1],it[3]] }
        .groupTuple(by: [0,1])
        .map{ it -> [it[0],it[1],it[2].flatten().collect()] }

    ch_input_filter_df_nv = MERGE_PROCESSED_VCF.out.df_nv
        .combine(ch_chosen_after_bulk_early, by: 0)

    FILTER_DF_NV_BY_CHOSEN_VARIANTS (
        ch_input_filter_df_nv,
        params.publish_dir,
        params.enable_publish
    )

    ch_input_create_tab_nvnr = FILTER_DF_NV_BY_CHOSEN_VARIANTS.out.df_nv_filtered
        .combine(ch_input_df_nr, by: 0)
        .map{ it -> [it[0], it[2], it[1], it[3].flatten().collect()] }

    CREATE_TAB_NVNR (
        ch_input_create_tab_nvnr,
        params.publish_dir,
        params.enable_publish
    )

    SEQUOIA_BINOM_BETABINOM_TAB_NV_NR (
        CREATE_TAB_NVNR.out,
        params.aggregated_min_mean_depth,
        params.aggregated_max_mean_depth,
        params.gender,
        params.publish_dir,
        params.disable_publish
    )

    ch_input_concat_filter_bb = SEQUOIA_BINOM_BETABINOM_TAB_NV_NR.out.df_filter
        .map{ it -> [it[0],it[2]] }
        .groupTuple(by:[0])
        .map{ it -> [it[0],it[1].flatten().collect()] }

    CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR (
        ch_input_concat_filter_bb,
        params.first_pass_binomial_cutoff,
        params.first_pass_betabinomial_cutoff,
        params.publish_dir,
        params.enable_publish
    )

    ch_input_vep = MERGE_PROCESSED_VCF.out.merged_vcf
        .combine(CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR.out.chosen_variants, by: 0)

    SUBSET_MERGED_VCF_CHOSEN_VARIANTS (
        ch_input_vep,
        params.publish_dir,
        params.enable_publish
    )

    ch_subset_by_chr = SUBSET_MERGED_VCF_CHOSEN_VARIANTS.out.subset_vcf
        .combine(Channel.of( params.chrs instanceof List ? params.chrs : params.chrs.tokenize(',') ).flatMap())

    SPLIT_SUBSET_VCF_BY_CHR (
        ch_subset_by_chr,
        params.publish_dir,
        params.enable_publish
    )

    if (params.vep_cache_dir) {
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

        ch_vep_by_group = SORT_INDEX_VEP_VCF.out.vep_vcf_chr
            .groupTuple(by: 0)
            .map { group, list_chr, list_vcf, list_tbi ->
                [group, [list_chr, list_vcf].transpose().sort{ a, b -> a[0] <=> b[0] }.collect{ it[1] }]
            }

        MERGE_VEP_VCF_BY_GROUP (
            ch_vep_by_group,
            params.publish_dir,
            params.enable_publish
        )

        FILTER_VEP_GERMLINE (
            MERGE_VEP_VCF_BY_GROUP.out.vep_vcf.map { group, vcf, tbi -> tuple(group, [vcf]) },
            params.vep_max_af_1kg,
            params.vep_filter_by_existing_variation,
            params.publish_dir,
            params.enable_publish
        )
        ch_chosen_variants = FILTER_VEP_GERMLINE.out.chosen_variants
    } else {
        ch_chosen_variants = CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR.out.chosen_variants
    }

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
        .groupTuple(by: [0,1])
        .map{ it -> [it[0],it[1],it[2].flatten().collect()] }
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
        params.enable_publish
    )

    ch_input_tabs_group_level = CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES.out.tabs
        .map{ it -> [it[0],it[2],it[3]] }
        .groupTuple(by:0)
        .map{ it -> [it[0],it[1].flatten().collect(),it[2].flatten().collect()] }

    CUSTOM_CREATE_GROUP_LEVEL_TAB_DFS (
        ch_input_tabs_group_level,
        params.publish_dir,
        params.enable_publish
    )

    SEQUOIA_SECOND_FILTER (
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

    // ── CUSTOM_VARIANT_FILTER_PROVENANCE ──────────────────────────────────────
    ch_pileup_per_group = CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES.out.pileup
        .map { group, chr, f -> tuple(group, f) }
        .groupTuple()

    ch_tab_nvnr_per_group = CREATE_TAB_NVNR.out
        .map { group, chr, mat_nv, mat_nr -> tuple(group, [mat_nv, mat_nr]) }
        .groupTuple(by: 0)
        .map { group, list -> tuple(group, list.flatten()) }

    ch_input_variant_provenance = GET_VARIANTS_FROM_MERGED_VCF.out.all_variants
        .combine(FILTER_CHOSEN_VARIANTS_BY_BULK.out.bulk_filter_provenance, by: 0)
        .combine(CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR.out.res_df,        by: 0)
        .combine(FILTER_VEP_GERMLINE.out.filter_provenance,                  by: 0)
        .combine(ch_pileup_per_group,                                        by: 0)
        .combine(SEQUOIA_SECOND_FILTER.out.df_filter,                        by: 0)
        .combine(ch_tab_nvnr_per_group,                                      by: 0)

    CUSTOM_VARIANT_FILTER_PROVENANCE (
        ch_input_variant_provenance,
        params.publish_dir,
        params.enable_publish
    )

    // ── Per-sample VCF annotation ─────────────────────────────────────────────
    LIST_SAMPLES_FROM_GROUP_VCF (
        MERGE_VEP_VCF_BY_GROUP.out.vep_vcf,
        params.publish_dir,
        params.enable_publish
    )

    ch_per_sample = LIST_SAMPLES_FROM_GROUP_VCF.out.sample_list
        .map { group, f -> f.readLines().collect { s -> tuple(group, s) } }
        .flatMap { it }

    ch_annotate_input = ch_per_sample
        .combine(MERGE_VEP_VCF_BY_GROUP.out.vep_vcf,                        by: 0)
        .combine(CUSTOM_VARIANT_FILTER_PROVENANCE.out.vcf_annotation_table, by: 0)
        .combine(CUSTOM_VARIANT_FILTER_PROVENANCE.out.pileup_focal,         by: 0)

    ANNOTATE_SAMPLE_VCF (
        ch_annotate_input,
        params.publish_dir,
        params.enable_publish
    )

    // ── Genotype tables ───────────────────────────────────────────────────────
    GENOTYPE_TABLE_FROM_ANNOTATED_VCF (
        ANNOTATE_SAMPLE_VCF.out.annotated_vcf,
        params.publish_dir,
        params.enable_publish
    )

    // ── NR/NV matrices ────────────────────────────────────────────────────────
    ch_annotated_grouped = ANNOTATE_SAMPLE_VCF.out.annotated_vcf
        .groupTuple(by: 0)

    CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF (
        ch_annotated_grouped,
        params.publish_dir,
        params.enable_publish
    )

    // ── Phylogenetic trees ────────────────────────────────────────────────────
    SEQUOIA_PHYLOGENY_SNV (
        CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF.out.nr_nv_matrices,
        params.gender,
        params.sequoia_vaf_absent,
        params.sequoia_vaf_present,
        params.sequoia_tree_mut_pval,
        params.sequoia_keep_ancestral,
        params.sequoia_split_trees,
        params.sequoia_genotype_conv_prob,
        params.sequoia_min_pval_for_true_somatic,
        params.sequoia_min_variant_reads_shared,
        params.sequoia_min_vaf_shared,
        params.sequoia_create_multi_tree,
        params.sequoia_mpboot_path,
        params.publish_dir,
        params.enable_publish
    )

    SEQUOIA_PHYLOGENY_INDEL (
        CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF.out.nr_nv_matrices,
        params.gender,
        params.sequoia_vaf_absent,
        params.sequoia_vaf_present,
        params.sequoia_tree_mut_pval,
        params.sequoia_keep_ancestral,
        params.sequoia_split_trees,
        params.sequoia_genotype_conv_prob,
        params.sequoia_min_pval_for_true_somatic,
        params.sequoia_min_variant_reads_shared,
        params.sequoia_min_vaf_shared,
        params.sequoia_create_multi_tree,
        params.sequoia_mpboot_path,
        params.publish_dir,
        params.enable_publish
    )

    SEQUOIA_PHYLOGENY_BOTH (
        CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF.out.nr_nv_matrices,
        params.gender,
        params.sequoia_vaf_absent,
        params.sequoia_vaf_present,
        params.sequoia_tree_mut_pval,
        params.sequoia_keep_ancestral,
        params.sequoia_split_trees,
        params.sequoia_genotype_conv_prob,
        params.sequoia_min_pval_for_true_somatic,
        params.sequoia_min_variant_reads_shared,
        params.sequoia_min_vaf_shared,
        params.sequoia_create_multi_tree,
        params.sequoia_mpboot_path,
        params.publish_dir,
        params.enable_publish
    )

    // ── Tree topology comparison ──────────────────────────────────────────────
    TREES_COMPARE_SIMILARITIES (
        SEQUOIA_PHYLOGENY_SNV.out.phylogeny_outputs
            .join( SEQUOIA_PHYLOGENY_INDEL.out.phylogeny_outputs )
            .join( SEQUOIA_PHYLOGENY_BOTH.out.phylogeny_outputs ),
        params.publish_dir,
        params.enable_publish
    )

    // ── Variant placement ─────────────────────────────────────────────────────
    CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF.out.nr_nv_matrices
        .map { group, nr_list, nv_list ->
            def nrs = nr_list instanceof List ? nr_list : [nr_list]
            def nvs = nv_list instanceof List ? nv_list : [nv_list]
            tuple(group,
                  nrs.find { it.name =~ /unfiltered/ },
                  nvs.find { it.name =~ /unfiltered/ })
        }
        .multiMap { group, nr, nv ->
            snv:         tuple(group, nr, nv)
            indel:       tuple(group, nr, nv)
            both:        tuple(group, nr, nv)
            postprocess: tuple(group, nr, nv)
        }
        .set { ch_unfiltered }

    def findPileupTree = { outputs, subdir ->
        def out_list = outputs instanceof List ? outputs : [outputs]
        def pileup_dir = out_list.find { it.isDirectory() && it.name == subdir }
        if (!pileup_dir) return file('/dev/null')
        def names = pileup_dir.list()
        if (!names) return file('/dev/null')
        def tf_name = names.find { it.endsWith('.treefile') }
        return tf_name ? pileup_dir.resolve(tf_name) : file('/dev/null')
    }

    ch_snv_placement_input = SEQUOIA_PHYLOGENY_SNV.out.phylogeny_outputs
        .map { group, outputs -> tuple(group, findPileupTree(outputs, "output_snv_pileup")) }
        .combine(ch_unfiltered.snv, by: 0)

    ch_indel_placement_input = SEQUOIA_PHYLOGENY_INDEL.out.phylogeny_outputs
        .map { group, outputs -> tuple(group, findPileupTree(outputs, "output_indel_pileup")) }
        .combine(ch_unfiltered.indel, by: 0)

    ch_both_placement_input = SEQUOIA_PHYLOGENY_BOTH.out.phylogeny_outputs
        .map { group, outputs -> tuple(group, findPileupTree(outputs, "output_both_pileup")) }
        .combine(ch_unfiltered.both, by: 0)

    SEQUOIA_VARIANT_PLACEMENT_SNV (
        ch_snv_placement_input,
        params.gender,
        params.sequoia_vaf_absent,
        params.sequoia_vaf_present,
        params.sequoia_tree_mut_pval,
        params.sequoia_keep_ancestral,
        params.sequoia_create_multi_tree,
        params.sequoia_genotype_conv_prob,
        params.sequoia_min_pval_for_true_somatic,
        params.sequoia_min_variant_reads_shared,
        params.sequoia_min_vaf_shared,
        params.publish_dir,
        params.enable_publish
    )

    SEQUOIA_VARIANT_PLACEMENT_INDEL (
        ch_indel_placement_input,
        params.gender,
        params.sequoia_vaf_absent,
        params.sequoia_vaf_present,
        params.sequoia_tree_mut_pval,
        params.sequoia_keep_ancestral,
        params.sequoia_create_multi_tree,
        params.sequoia_genotype_conv_prob,
        params.sequoia_min_pval_for_true_somatic,
        params.sequoia_min_variant_reads_shared,
        params.sequoia_min_vaf_shared,
        params.publish_dir,
        params.enable_publish
    )

    SEQUOIA_VARIANT_PLACEMENT_BOTH (
        ch_both_placement_input,
        params.gender,
        params.sequoia_vaf_absent,
        params.sequoia_vaf_present,
        params.sequoia_tree_mut_pval,
        params.sequoia_keep_ancestral,
        params.sequoia_create_multi_tree,
        params.sequoia_genotype_conv_prob,
        params.sequoia_min_pval_for_true_somatic,
        params.sequoia_min_variant_reads_shared,
        params.sequoia_min_vaf_shared,
        params.publish_dir,
        params.enable_publish
    )

    // ── VAF + digital heatmaps ────────────────────────────────────────────────
    ch_gt_per_group = GENOTYPE_TABLE_FROM_ANNOTATED_VCF.out.genotype_table
        .map  { group, sample_name, tsv -> tuple(group, tsv) }
        .groupTuple(by: 0)

    ch_postprocess_input = SEQUOIA_VARIANT_PLACEMENT_SNV.out.placement_outputs
        .join(SEQUOIA_VARIANT_PLACEMENT_INDEL.out.placement_outputs, by: 0)
        .join(SEQUOIA_VARIANT_PLACEMENT_BOTH.out.placement_outputs,  by: 0)
        .combine(ch_unfiltered.postprocess, by: 0)
        .combine(ch_gt_per_group, by: 0)

    POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE (
        ch_postprocess_input,
        params.publish_dir,
        params.enable_publish
    )

    // ── Mutational signature analysis (optional) ──────────────────────────────
    ch_sig_zero_png   = Channel.value(file('/dev/null'))
    ch_sig_cosine_png = Channel.value(file('/dev/null'))

    if ( params.run_sigprofiler_assignment ) {

        SUBSET_ANNOTATED_VCFS_FOR_MUTSIG (
            ANNOTATE_SAMPLE_VCF.out.annotated_vcf,
            params.publish_dir,
            params.enable_publish
        )

        SIGPROFILER_ASSIGNMENT (
            SUBSET_ANNOTATED_VCFS_FOR_MUTSIG.out.filtered_vcf,
            params.sig_genome_build,
            params.sig_context_type,
            params.sig_exome,
            params.sig_export_probabilities,
            params.sig_export_probabilities_per_mutation,
            params.publish_dir,
            params.enable_publish
        )

        ch_all_assignment_dirs = SIGPROFILER_ASSIGNMENT.out.assignment_output
            .map { sample_name, dir -> dir }
            .collect()

        MERGE_SIGNATURE_ACTIVITIES (
            ch_all_assignment_dirs,
            params.publish_dir,
            params.enable_publish
        )

        if ( params.sig_filter_enabled ) {

            ch_filter_input = SUBSET_ANNOTATED_VCFS_FOR_MUTSIG.out.filtered_vcf
                .join( SIGPROFILER_ASSIGNMENT.out.assignment_output, by: 0 )

            FILTER_VARIANTS_BY_SIG_PROBABILITY (
                ch_filter_input,
                params.reference,
                params.sig_genome_build,
                params.sig_cosmic_version,
                params.publish_dir,
                params.enable_publish
            )

            ch_all_cosine_summaries = FILTER_VARIANTS_BY_SIG_PROBABILITY.out.cosine_summary
                .collect()

            MERGE_SIG_COSINE_SIMILARITIES (
                ch_all_cosine_summaries,
                params.publish_dir,
                params.enable_publish
            )

            FILTER_ACTIVITIES_BY_COSINE (
                MERGE_SIG_COSINE_SIMILARITIES.out.cosine_matrix,
                MERGE_SIGNATURE_ACTIVITIES.out.merged_activities,
                params.sig_min_cosine_for_activity,
                params.publish_dir,
                params.enable_publish
            )

            PLOT_ZERO_FILTERED_SIGNATURE_ACTIVITIES (
                FILTER_ACTIVITIES_BY_COSINE.out.zero_filtered_activities,
                "zero_filtered",
                params.publish_dir,
                params.enable_publish
            )

            PLOT_COSINE_FILTERED_SIGNATURE_ACTIVITIES (
                FILTER_ACTIVITIES_BY_COSINE.out.filtered_activities,
                "cosine_filtered",
                params.publish_dir,
                params.enable_publish
            )

            ch_sig_zero_png   = PLOT_ZERO_FILTERED_SIGNATURE_ACTIVITIES.out.signature_plot
            ch_sig_cosine_png = PLOT_COSINE_FILTERED_SIGNATURE_ACTIVITIES.out.signature_plot
        }
    }

    // ── Germline identification and per-sample subsetting ────────────────────
    ch_identify_germline_input = MERGE_PROCESSED_VCF.out.merged_vcf
        .combine(CUSTOM_VARIANT_FILTER_PROVENANCE.out.master_table,         by: 0)
        .combine(FILTER_VEP_GERMLINE.out.filter_provenance,                  by: 0)
        .combine(FILTER_CHOSEN_VARIANTS_BY_BULK.out.bulk_filter_provenance,  by: 0)

    IDENTIFY_GERMLINE_FROM_STATS (
        ch_identify_germline_input,
        params.germline_prev_pct,
        params.publish_dir,
        params.enable_publish
    )

    EXTRACT_GERMLINE_PREVALENCE_TABLE (
        IDENTIFY_GERMLINE_FROM_STATS.out.annotated_vcf,
        params.publish_dir,
        params.enable_publish
    )

    PLOT_GERMLINE_PREVALENCE_DISTRIBUTIONS (
        EXTRACT_GERMLINE_PREVALENCE_TABLE.out.prevalence_table,
        params.germline_prev_pct,
        params.publish_dir,
        params.enable_publish
    )

    ch_subset_germline_input = ch_per_sample
        .combine(IDENTIFY_GERMLINE_FROM_STATS.out.annotated_vcf, by: 0)

    SUBSET_MERGED_VCF_HIGH_CONFIDENCE_GERMLINE_FROM_STATS (
        ch_subset_germline_input,
        params.publish_dir,
        params.enable_publish
    )

    // ── ADO analysis on per-sample germline VCFs ──────────────────────────────
    ch_ado_input = SUBSET_MERGED_VCF_HIGH_CONFIDENCE_GERMLINE_FROM_STATS.out.stats_vcf
        .map { group, sample_name, vcf, tbi -> tuple("${sample_name}_stats", vcf, tbi) }
        .mix(
            SUBSET_MERGED_VCF_HIGH_CONFIDENCE_GERMLINE_FROM_STATS.out.vep_vcf
                .map { group, sample_name, vcf, tbi -> tuple("${sample_name}_vep", vcf, tbi) },
            SUBSET_MERGED_VCF_HIGH_CONFIDENCE_GERMLINE_FROM_STATS.out.bulk_vcf
                .map { group, sample_name, vcf, tbi -> tuple("${sample_name}_bulk", vcf, tbi) }
        )

    CREATE_ADO_TABLE_FROM_GERMLINE_VCF (
        ch_ado_input,
        params.ado_sample_prop,
        params.publish_dir,
        params.enable_publish
    )

    SUMMARIZE_ADO_INTERVALS (
        CREATE_ADO_TABLE_FROM_GERMLINE_VCF.out.ado_table
            .filter { sample_name, tsv -> tsv.size() > 0 },
        params.ado_cov_cutoff,
        params.publish_dir,
        params.enable_publish
    )

    ch_df_sum_by_prov = SUMMARIZE_ADO_INTERVALS.out.df_sum
        .flatten()
        .branch {
            stats: it.name.contains('_stats')
            vep:   it.name.contains('_vep')
            bulk:  it.name.contains('_bulk')
        }

    CONCAT_ADO_STATS ( ch_df_sum_by_prov.stats.collect(), 'stats', params.publish_dir, params.enable_publish )
    CONCAT_ADO_VEP   ( ch_df_sum_by_prov.vep.collect(),   'vep',   params.publish_dir, params.enable_publish )
    CONCAT_ADO_BULK  ( ch_df_sum_by_prov.bulk.collect(),  'bulk',  params.publish_dir, params.enable_publish )

    ch_plot_ado_tables = CONCAT_ADO_STATS.out.merged_ADO
        .mix( CONCAT_ADO_VEP.out.merged_ADO  )
        .mix( CONCAT_ADO_BULK.out.merged_ADO )
        .collect()

    ch_plot_ado_summaries = CONCAT_ADO_STATS.out.summary_ADO
        .mix( CONCAT_ADO_VEP.out.summary_ADO  )
        .mix( CONCAT_ADO_BULK.out.summary_ADO )
        .collect()

    PLOT_ADO_GERMLINE_COMPARISON (
        ch_plot_ado_tables,
        ch_plot_ado_summaries,
        params.publish_dir,
        params.enable_publish
    )

    // ── Master report ─────────────────────────────────────────────────────────
    ch_postprocess_grouped = POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE.out.postprocess_outputs
        .map { group, snv_dir, indel_dir, both_dir -> tuple(group, [snv_dir, indel_dir, both_dir]) }

    ch_compile_input = ch_postprocess_grouped
        .join(
            TREES_COMPARE_SIMILARITIES.out.comparison_results
                .map { group, dirs, pdf -> tuple(group, pdf) },
            by: 0
        )
        .join( CUSTOM_VARIANT_FILTER_PROVENANCE.out.combined_report, by: 0 )
        .join( PLOT_GERMLINE_PREVALENCE_DISTRIBUTIONS.out.prevalence_plot, by: 0 )

    COMPILE_MASTER_REPORT (
        ch_compile_input,
        ch_sig_zero_png,
        ch_sig_cosine_png,
        PLOT_ADO_GERMLINE_COMPARISON.out.combined_plot,
        params.publish_dir,
        params.enable_publish
    )

}
