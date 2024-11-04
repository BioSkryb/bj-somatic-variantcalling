nextflow.enable.dsl=2

// IMPORT MODULES
include { SENTIEON_DRIVER_TNSEQ } from '../../modules/sentieon/driver/tnseq/main.nf' addParams( timestamp: params.timestamp )
include { PSEUDO_BULK_WF } from '../../subworkflows/pseudobulk_sc_wf/main.nf' addParams( timestamp: params.timestamp )
include { SENTIEON_DRIVER_TNSCOPE; SENTIEON_DRIVER_TNSCOPE_PANEL_OF_NORMAL; SENTIEON_DRIVER_TNSCOPE_TUMOR_PANEL_OF_NORMAL; GENERATE_PANEL_OF_NORMAL } from '../../modules/sentieon/driver/tnscope/main.nf' addParams( timestamp: params.timestamp )

params.reference                    = params.genomes [ params.genome ] [ 'reference' ]
params.dbsnp                        = params.genomes [ params.genome ] [ 'dbsnp' ]
params.dbsnp_index                  = params.genomes [ params.genome ] [ 'dbsnp_index' ]


workflow SOMATIC_VARIANT_WORKFLOW_MATCHNORMAL {
    take:
        ch_input
        ch_reference
        ch_dbsnp
        ch_dbsnp_index
        ch_interval
        ch_somatic_variant_caller
        ch_publish_dir
        ch_enable_publish

    main:
        if (ch_somatic_variant_caller == 'tnscope') {
            SENTIEON_DRIVER_TNSCOPE ( 
                                        ch_input,
                                        ch_reference,
                                        ch_dbsnp,
                                        ch_dbsnp_index,
                                        ch_interval,
                                        ch_publish_dir,
                                        ch_enable_publish
                                        )
            emit:
                vcf = SENTIEON_DRIVER_TNSCOPE.out.vcf
                version = SENTIEON_DRIVER_TNSCOPE.out.version
                vcf_only = SENTIEON_DRIVER_TNSCOPE.out.vcf_only
                vcf_index_only = SENTIEON_DRIVER_TNSCOPE.out.vcf_index_only
        } else {
            SENTIEON_DRIVER_TNSEQ ( 
                                    ch_input,
                                    ch_reference,
                                    ch_publish_dir,
                                    ch_enable_publish
                                    )
            emit:
                vcf = SENTIEON_DRIVER_TNSEQ.out.vcf
                version = SENTIEON_DRIVER_TNSEQ.out.version
        
        }
    emit:
        vcf = vcf
        version = version

}

workflow SOMATIC_VARIANT_WORKFLOW_PSEUDOBULK {
    take:

        pseudobulk_input_bam
        ch_genome
        ch_reference
        ch_dbsnp
        ch_dbsnp_index
        ch_mills
        ch_mills_index
        ch_onekg
        ch_onekg_index
        ch_reads
        ch_pseudobulk_reads
        ch_bam_recaltab
        ch_interval
        ch_somatic_variant_caller
        ch_publish_dir
        ch_enable_publish

    main:

        PSEUDO_BULK_WF (
                            pseudobulk_input_bam,
                            ch_genome,
                            ch_reference,
                            ch_dbsnp,
                            ch_dbsnp_index,
                            ch_mills,
                            ch_mills_index,
                            ch_onekg,
                            ch_onekg_index,
                            params.publish_dir,
                            params.enable_publish
                        )

        // Rearranging the meta index so that the pseudo sample name (which is group name) comes in line with group in the ch_pseudobulk_reads
        pseudo_bam_recal_rearranged = PSEUDO_BULK_WF.out.pseudo_bam_recal
                                        .map { sample, bam, bai, recal -> [bam, bai, sample, recal]}

        ch_bam_bulk = pseudo_bam_recal_rearranged.join(ch_pseudobulk_reads, by: 2)
            .map { pseudo_group, bam, bai, recal_table, sample_name, reads, isbulk, n_reads -> [pseudo_group, [bam, bai], recal_table, pseudo_group] }

        // Setting the single/tumor channel using the samples which has isbulk set as false
        ch_reads_single = ch_reads.filter { !it[3].toBoolean() }
        ch_bam_single = ch_bam_recaltab.join(ch_reads_single).map 
            { sample_name, bam, bai, recal_table, fasta, groups, isbulk -> [sample_name, [bam, bai], recal_table, groups] }

        ch_input = ch_bam_single.combine( ch_bam_bulk, by: 3 )

        if (ch_somatic_variant_caller == 'tnscope') {
            SENTIEON_DRIVER_TNSCOPE ( 
                                        ch_input,
                                        ch_reference,
                                        ch_dbsnp,
                                        ch_dbsnp_index,
                                        ch_interval,
                                        ch_publish_dir,
                                        ch_enable_publish
                                        )
            emit:
                vcf = SENTIEON_DRIVER_TNSCOPE.out.vcf
                version = SENTIEON_DRIVER_TNSCOPE.out.version
                vcf_only = SENTIEON_DRIVER_TNSCOPE.out.vcf_only
                vcf_index_only = SENTIEON_DRIVER_TNSCOPE.out.vcf_index_only
        } else {
            SENTIEON_DRIVER_TNSEQ ( 
                                    ch_input,
                                    ch_reference,
                                    ch_publish_dir,
                                    ch_enable_publish
                                    )
            emit:
                vcf = SENTIEON_DRIVER_TNSEQ.out.vcf
                version = SENTIEON_DRIVER_TNSEQ.out.version
        
        }
    emit:
        vcf = vcf
        version = version
        picard_version = PSEUDO_BULK_WF.out.picard_version
        samtools_version = PSEUDO_BULK_WF.out.samtools_version

}

workflow SOMATIC_VARIANT_WORKFLOW_PANELNORMAL {
    take:
        ch_panel_of_normal_bam
        ch_single_bam
        ch_reference
        ch_dbsnp
        ch_dbsnp_index
        ch_interval
        ch_panel_of_normal_vcf
        ch_publish_dir
        ch_enable_publish

    main:
        
        // Create panel of normal if not exist
        if ( ch_panel_of_normal_vcf == "") {
            ch_panel_of_normal_bam.view()

            SENTIEON_DRIVER_TNSCOPE_PANEL_OF_NORMAL ( ch_panel_of_normal_bam, 
                                                    ch_reference, 
                                                    ch_interval, 
                                                    ch_publish_dir, 
                                                    ch_enable_publish ) 

            GENERATE_PANEL_OF_NORMAL ( SENTIEON_DRIVER_TNSCOPE_PANEL_OF_NORMAL.out.vcf_only.collect(),
                                       SENTIEON_DRIVER_TNSCOPE_PANEL_OF_NORMAL.out.vcf_index_only.collect() )

            // combine the tumor to panel of normal pairs
            ch_tumor = ch_single_bam.combine(GENERATE_PANEL_OF_NORMAL.out)
        
        } else {
            panel_of_normal_ch = Channel.fromPath(ch_panel_of_normal_vcf)
                                    .map { it -> [ it, it + ".tbi" ] }
            panel_of_normal_ch.view()
            ch_tumor = ch_single_bam.combine(panel_of_normal_ch)
        }

        ch_tumor.view()

        // run tnscope with generic panel of normal file
        SENTIEON_DRIVER_TNSCOPE_TUMOR_PANEL_OF_NORMAL ( ch_tumor, 
                                    ch_reference, 
                                    ch_dbsnp, 
                                    ch_dbsnp_index, 
                                    ch_interval, 
                                    ch_publish_dir, 
                                    ch_enable_publish
                                    )
            
    emit:
        vcf = SENTIEON_DRIVER_TNSCOPE_TUMOR_PANEL_OF_NORMAL.out.vcf
        version =  SENTIEON_DRIVER_TNSCOPE_TUMOR_PANEL_OF_NORMAL.out.version

}
