nextflow.enable.dsl=2

// IMPORT MODULES

include { CUSTOM_BAM_SUBSAMPLE } from '../../modules/bioskryb/custom_bam_subsample/main.nf'
include { CUSTOM_BAM_CONCATENATE_FILES } from '../../modules/bioskryb/custom_bam_concatenate_files/main.nf'
include { PICARD_ADDORREPLACEREADGROUPS } from '../../modules/picard/addorreplacereadgroups/main.nf'
include { SENTIEON_ALGORITHM } from '../../modules/sentieon/driver/alignment/main.nf'

workflow PSEUDO_BULK_WF {

    take:
        ch_bam
        ch_genome
        ch_reference
        ch_dbsnp
        ch_dbsnp_index
        ch_mills
        ch_mills_index
        ch_onekg
        ch_onekg_index
        ch_publish_dir
        ch_enable_publish


    main:

        ch_bam_subsample = ch_bam.map { sample, bam, bai, groups, n_reads -> [sample, bam, bai, n_reads] }

        ch_groups = ch_bam.map { sample, bam, bai, groups, n_reads -> [sample,groups] }

        CUSTOM_BAM_SUBSAMPLE(
                                ch_bam_subsample,
                                params.publish_dir,
                                params.enable_publish
                            )

        ch_sample_with_groups = ch_groups.join(CUSTOM_BAM_SUBSAMPLE.out.bam_with_samplename)
        ch_sample_with_groups.groupTuple(by: 1)
                .map { samples, group, bams, bais ->
                    return [group, (bams + bais)]
                }
                .set { bams_to_concat }

        bams_to_concat.toList().collect { list ->
                            list.collect { group, files -> [group, files.size() / 2] }}
                            .subscribe { sizes ->
                                    sizes.each { group, size ->
                                        if (size < 3) {
                                            throw new RuntimeException("ERROR: The number of BAMs in group ${group} is less than 3 (${size}), which is insufficient to create pseudobulk. Exiting.")
                                        } else {
                                            println "Number of BAMs in group ${group} available for pseudobulk creation: ${size}"
                                        }
                                    }
                                } 

        CUSTOM_BAM_CONCATENATE_FILES (

                                    bams_to_concat,
                                    params.publish_dir,
                                    params.disable_publish
        )


        PICARD_ADDORREPLACEREADGROUPS (

                                    CUSTOM_BAM_CONCATENATE_FILES.out.bam,
                                    params.publish_dir,
                                    params.disable_publish

        )

        SENTIEON_ALGORITHM (
                                ch_genome,
                                PICARD_ADDORREPLACEREADGROUPS.out.bam,
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

    emit:
        pseudo_bam = PICARD_ADDORREPLACEREADGROUPS.out.bam
        pseudo_bam_recal = SENTIEON_ALGORITHM.out.bam_recal_table
        picard_version = PICARD_ADDORREPLACEREADGROUPS.out.version
        samtools_version = CUSTOM_BAM_SUBSAMPLE.out.version

}