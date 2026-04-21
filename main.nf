nextflow.enable.dsl=2

include { SOMATIC_SNP_INDEL_FILTERING_WF } from './nf-bioskryb-utils/subworkflows/somatic_heuristic_filter/main.nf'

params.reference = params.genomes [ params.genome ] [ 'reference' ]

workflow {

    if ( !params.input_csv ) {
        error "ERROR: input_csv is not defined.\nPlease add --input_csv <path/to/input.csv> to specify the input CSV file."
    }

    if ( !params.publish_dir ) {
        error "ERROR: publish_dir is not defined.\nPlease add --publish_dir s3://<bucket_name>/<project_name> to specify where the pipeline outputs will be stored."
    }

    SOMATIC_SNP_INDEL_FILTERING_WF()

}
