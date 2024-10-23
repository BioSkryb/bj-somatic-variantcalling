nextflow.enable.dsl=2

// IMPORT MODULES
include { SENTIEON_DRIVER_TNSCOPE } from '../../modules/sentieon/driver/tnscope/main.nf' addParams( timestamp: params.timestamp )
include { SENTIEON_DRIVER_TNSEQ } from '../../modules/sentieon/driver/tnseq/main.nf' addParams( timestamp: params.timestamp )

params.reference                    = params.genomes [ params.genome ] [ 'reference' ]
params.dbsnp                        = params.genomes [ params.genome ] [ 'dbsnp' ]
params.dbsnp_index                  = params.genomes [ params.genome ] [ 'dbsnp_index' ]


workflow SOMATIC_VARIANT_WORKFLOW {
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
