nextflow.enable.dsl=2
import groovy.json.JsonOutput

include { printHeader; helpMessage } from './help' params ( params )
include { PUBLISH_INPUT_DATASET_WF } from './nf-bioskryb-utils/modules/bioskryb/publish_input_dataset/main.nf' addParams(timestamp: params.timestamp)
include { SOMATIC_VARIANT_WORKFLOW } from './nf-bioskryb-utils/subworkflows/somatic_variant_calling/main.nf' params ( params )
include { VARIANT_ANNOTATION_WF } from './nf-bioskryb-utils/subworkflows/variant_annotation/main.nf' params ( params )
include { SENTIEON_DRIVER_TNSCOPE } from './nf-bioskryb-utils/modules/sentieon/driver/tnscope/main.nf' addParams( timestamp: params.timestamp )
include { SENTIEON_DRIVER_TNSEQ } from './nf-bioskryb-utils/modules/sentieon/driver/tnseq/main.nf' addParams( timestamp: params.timestamp )
include { SENTIEON_DRIVER_METRICS_WF } from './nf-bioskryb-utils/modules/sentieon/driver/metrics/main.nf' addParams(timestamp: params.timestamp)
include { SENTIEON_DRIVER_COVERAGEMETRICS_WF } from './nf-bioskryb-utils/modules/sentieon/driver/coveragemetrics/main.nf' addParams(timestamp: params.timestamp)
include { CUSTOM_DATA_PROCESSING_WF } from './modules/local/custom_data_processing/main.nf' addParams(timestamp: params.timestamp)
include { CUSTOM_REPORT_WF } from './modules/local/custom_report/main.nf' addParams(timestamp: params.timestamp)
include { MULTIQC_WF } from './nf-bioskryb-utils/modules/multiqc/main.nf' addParams(timestamp: params.timestamp)
include { SET_PSEUDO_BULK } from './modules/local/set_pseudo_bulk/main.nf' addParams(timestamp: params.timestamp)
include { CUSTOM_READ_COUNTS_WF } from './nf-bioskryb-utils/modules/bioskryb/custom_read_counts/main.nf' addParams(timestamp: params.timestamp)
include { SENTIEON_ALIGNMENT } from './nf-bioskryb-utils/modules/sentieon/driver/alignment/main.nf' addParams( timestamp: params.timestamp )
include { PSEUDO_BULK_WF } from './nf-bioskryb-utils/subworkflows/pseudobulk_sc_wf/main.nf' addParams( timestamp: params.timestamp )
include { REPORT_VERSIONS_WF } from './nf-bioskryb-utils/modules/bioskryb/report_tool_versions/main.nf' addParams(timestamp: params.timestamp)


if ( params.help ) {
    helpMessage()
    exit 0
}


workflow {
    
    printHeader()
    
    // Check if publishDir specified
    
    if ( params.publish_dir == "") {
        exit 1, "ERROR: publish_dir is not defined.\nPlease add --publish_dir s3://<bucket_name>/<project_name> to specify where the pipeline outputs will be stored."
    }
    ch_input_csv = Channel.fromPath(params.input_csv, checkIfExists: true).collect()

    if ( !params.is_bam ){
        ch_reads = Channel.fromPath( params.input_csv ).splitCsv( header:true )
                .map { row -> [ row.biosampleName, [ row.read1, row.read2 ], row.groups, row.isbulk ] }
    } else {
        ch_reads = Channel.fromPath( params.input_csv ).splitCsv( header:true )
                .map { row -> [ row.biosampleName, [ row.bam, row.bam + ".bai", row.bam.replace(".bam", "_recal_data.table"), row.bam.replace(".bam", ".dedup_sentieonmetrics.txt") ], row.groups, row.isbulk ] }
    }
        
    ch_reads.ifEmpty{ exit 1, "ERROR: Input csv file is empty." }

    if ( params.mode == 'wgs' ) {
            params.calling_intervals_filename   = params.genomes [ params.genome ] [ 'calling_intervals_filename' ]
    }

    if ( params.mode == 'exome' ) {
        params.calling_intervals_filename   = params.genomes [ params.genome ] [ params.exome_panel ] [ 'calling_intervals_filename' ]
    }

    /*
    ========================================================================================
        PROCESS PRIMARY DATA
    ========================================================================================
    */
    PUBLISH_INPUT_DATASET_WF (
                                        ch_input_csv,
                                        params.publish_dir,
                                        params.enable_publish
                                    )

    /*
    ========================================================================================
        SET_PSEUDO_BULK
    ========================================================================================
    */

    // Counting the reads for each sample
    CUSTOM_READ_COUNTS_WF (
                            ch_reads.map{ [it[0], it[1]] },
                            params.publish_dir,
                            params.enable_publish
                        )
    // Process to check if bulk samples are present and choose samples for pseudobulk generation    
    SET_PSEUDO_BULK(ch_input_csv, CUSTOM_READ_COUNTS_WF.out.combined_read_counts, params.publish_dir, params.enable_publish)
    SET_PSEUDO_BULK.out.bulk_bool.map { it.text.trim() }.set { bulk_bool_value }

    if ( !params.is_bam ){
        ch_pseudobulk_reads = SET_PSEUDO_BULK.out.pseudo_bulk_samples.splitCsv( header:true )
                    .map { row -> [ row.biosampleName, [ row.read1, row.read2 ], row.groups, row.isbulk, row.n_reads ] }
    } else {  
        ch_pseudobulk_reads = SET_PSEUDO_BULK.out.pseudo_bulk_samples.splitCsv( header:true )
                    .map { row -> [ row.biosampleName, [ row.bam, row.bam + ".bai", row.bam.replace(".bam", "_recal_data.table")], row.groups, row.isbulk, row.n_reads ] }

        ch_bam_recaltab = ch_reads.map{ sample, files, groups, isbulk -> [sample, files[0], files[1], files[2]] }
        ch_bam_only = ch_reads.map{ sample, files, groups, isbulk -> [sample, files[0], files[1]] }
        ch_dedup_metrics = ch_reads.map{ sample, files, groups, isbulk -> [sample, files[3]] }
    }
    ch_pseudobulk_reads_list = ch_pseudobulk_reads.toList()
    ch_pseudobulk_reads_list.view()
        

    // Checks before pseudobulk generation
    bulk_bool_value.subscribe { bulk_bool ->
        if (!bulk_bool.toBoolean()) {
            println("No bulk samples found. Proceeding with pseudobulk generation.")
            ch_pseudobulk_reads_list.subscribe { reads ->
                if (!params.pseudo_bulk_run || "${params.pseudo_bulk_run}".equalsIgnoreCase("false")){
                    log.error("ERROR: The parameter 'pseudo_bulk_run' must be set to 'true' when there are no bulk samples in the input.csv. Exiting.")
                    exit 1
                }
            }
        }
    }


    if ( !params.is_bam ){
        /*
        ========================================================================================
            ALIGNMENT
        ========================================================================================
        */
        SENTIEON_ALIGNMENT ( 
                                params.genome,
                                ch_reads.map{ [it[0], it[1]] },
                                params.reference,
                                params.dbsnp,
                                params.dbsnp_index,
                                params.mills,
                                params.mills_index,
                                params.onekg,
                                params.onekg_index,
                                params.is_fasterq,
                                params.publish_dir,
                                params.enable_publish
                            )
        
        pseudobulk_input_bam = SENTIEON_ALIGNMENT.out.bam.join(ch_pseudobulk_reads, by: 0)
        .map { sample, bam, bai, reads, groups, isbulk, n_reads -> [sample, bam, bai, groups, n_reads] }

        ch_bam_recaltab = SENTIEON_ALIGNMENT.out.bam_recal_table
        ch_bam_only = SENTIEON_ALIGNMENT.out.bam
        ch_dedup_metrics = SENTIEON_ALIGNMENT.out.dedup_metrics

        // Collecting the bam files for the output
        SENTIEON_ALIGNMENT.out.bam
                .collectFile( name: "bam_files.txt", newLine: true, sort: { it[0] }, storeDir: "${params.tmp_dir}" )
                    { it[0] + "\t" + "${params.publish_dir}_${params.timestamp}/secondary_analyses/alignment/" + it[1].getName() }

    } else{
        pseudobulk_input_bam = ch_bam_only.join(ch_pseudobulk_reads, by: 0)
        .map { sample, bam, bai, files, groups, isbulk, n_reads -> [sample, bam, bai, groups, n_reads] }

        ch_reads.map{ sample, files, groups, isbulk -> [sample, files[0], files[1]] }
        .collectFile( name: "bam_files.txt", newLine: true, sort: { it[0] }, storeDir: "${params.tmp_dir}" )
            { it[0] + "\t" + "${params.publish_dir}_${params.timestamp}/secondary_analyses/alignment/" + new File(it[1]).getName() }

    }
    

    /*
    ========================================================================================
        PSEUDOBULK
    ========================================================================================
    */

    if (params.pseudo_bulk_run || "${params.pseudo_bulk_run}".equalsIgnoreCase("true")){
        // Pseudobulk run and setting the channels for somatic variant calling
        PSEUDO_BULK_WF (
                            pseudobulk_input_bam,
                            params.genome,
                            params.reference,
                            params.dbsnp,
                            params.dbsnp_index,
                            params.mills,
                            params.mills_index,
                            params.onekg,
                            params.onekg_index,
                            params.publish_dir,
                            params.enable_publish
                        )

        ch_picard_version = PSEUDO_BULK_WF.out.picard_version
        ch_samtools_version = PSEUDO_BULK_WF.out.samtools_version

        // Rearranging the meta index so that the pseudo sample name (which is group name) comes in line with group in the ch_pseudobulk_reads
        pseudo_bam_recal_rearranged = PSEUDO_BULK_WF.out.pseudo_bam_recal
                                        .map { sample, bam, bai, recal -> [bam, bai, sample, recal]}

        ch_bam_bulk = pseudo_bam_recal_rearranged.join(ch_pseudobulk_reads, by: 2)
            .map { pseudo_group, bam, bai, recal_table, sample_name, reads, isbulk, n_reads -> [pseudo_group, [bam, bai], recal_table, pseudo_group] }
        
    } else{
        // Setting the bulk channel using the isbulk set to true samples
        ch_reads_bulk = ch_reads.filter { it[3].toBoolean() }
        ch_bam_bulk = ch_bam_recaltab.join(ch_reads_bulk).map 
            { sample_name, bam, bai, recal_table, fasta, groups, isbulk -> [sample_name, [bam, bai], recal_table, groups] }

        ch_picard_version = Channel.empty()
        ch_samtools_version = Channel.empty()

    }

    // Setting the single/tumor channel using the samples which has isbulk set as false
    ch_reads_single = ch_reads.filter { !it[3].toBoolean() }
    ch_bam_single = ch_bam_recaltab.join(ch_reads_single).map 
        { sample_name, bam, bai, recal_table, fasta, groups, isbulk -> [sample_name, [bam, bai], recal_table, groups] }

    /*
    ========================================================================================
        SOMATIC VARIANT WORKFLOW
    ========================================================================================
    */

    ch_input = ch_bam_single.combine( ch_bam_bulk, by: 3 )

    SOMATIC_VARIANT_WORKFLOW (
            ch_input,
            params.reference,
            params.dbsnp,
            params.dbsnp_index,
            params.calling_intervals_filename,
            params.somatic_variant_caller,
            params.publish_dir,
            params.enable_publish
        )

    ch_vep_version = Channel.empty()
    ch_bcftools_version = Channel.empty()
    ch_gatk_version = Channel.empty()
    ch_variantqc_version = Channel.empty()

    if ( !params.skip_variant_annotation ){

        /*
        ========================================================================================
            VARIANT ANNOTATION WORKFLOW
        ========================================================================================
        */

        VARIANT_ANNOTATION_WF (
                SOMATIC_VARIANT_WORKFLOW.out.vcf,
                params.genome,
                params.reference,
                params.vep_cache,
                params.publish_dir,
                params.enable_publish,
                params.disable_publish
            )

        ch_vep_version = VARIANT_ANNOTATION_WF.out.vep_version
        ch_bcftools_version = VARIANT_ANNOTATION_WF.out.bcftools_version
        ch_gatk_version = VARIANT_ANNOTATION_WF.out.gatk_version
        ch_variantqc_version = VARIANT_ANNOTATION_WF.out.variantqc_version

    }


    /*
    ========================================================================================
        METRICS
    ========================================================================================
    */

    SENTIEON_DRIVER_METRICS_WF (
                                    ch_bam_recaltab,
                                    params.reference,
                                    params.base_metrics_intervals,
                                    params.wgs_or_target_intervals,
                                    params.mode,
                                    params.publish_dir,
                                    params.enable_publish
                                )
    
    ch_sentieon_metrics_version = SENTIEON_DRIVER_METRICS_WF.out.version
    ch_sentieon_metrics = SENTIEON_DRIVER_METRICS_WF.out.metrics
                        
    if ( params.genome == "GRCh38" ){

        if ( !params.skip_gene_coverage ) {
    
            SENTIEON_DRIVER_COVERAGEMETRICS_WF (
                                                ch_bam_recaltab,
                                                params.reference,
                                                params.roi_intervals,
                                                params.refseq,
                                                params.publish_dir,
                                                params.enable_publish
                                        )
        }
    }

    ch_custom_data_processing_input = ch_dedup_metrics
                                                .combine (
                                                            ch_bam_only
                                                                .join( SENTIEON_DRIVER_METRICS_WF.out.metrics_tuple ),
                                                            by: 0
                                                         )


    CUSTOM_DATA_PROCESSING_WF ( 
                                ch_custom_data_processing_input,
                                params.mode,
                                params.publish_dir,
                                params.disable_publish
                                )

    ch_merge_metrics_report = CUSTOM_DATA_PROCESSING_WF.out.metrics
    ch_merge_metrics_version = CUSTOM_DATA_PROCESSING_WF.out.version
    ch_ado_summary_table = Channel.fromPath("$projectDir/assets/ado_dummy_file.txt", checkIfExists: true).collect()
    


    CUSTOM_REPORT_WF ( 
                        ch_merge_metrics_report.collect(),
                        ch_input_csv,
                        ch_ado_summary_table,
                        params.publish_dir,
                        params.enable_publish
                        )

    ch_custom_report = CUSTOM_REPORT_WF.out.mqc
    ch_custom_report_version = CUSTOM_REPORT_WF.out.version

    /*
    ========================================================================================
        REPORTING
    ========================================================================================
    */

    ch_tool_versions = ch_sentieon_metrics_version.take(1)
                                    .combine(ch_picard_version.take(1).ifEmpty([]))
                                    .combine(ch_samtools_version.take(1).ifEmpty([]))
                                    .combine((ch_vep_version ?: Channel.empty()).ifEmpty([]).first())
                                    .combine((ch_gatk_version ?: Channel.empty()).ifEmpty([]).first())
                                    .combine((ch_bcftools_version ?: Channel.empty()).ifEmpty([]).first())
                                    .combine((ch_variantqc_version ?: Channel.empty()).ifEmpty([]).first())
                                    
                                    
                                                
                                                    
    REPORT_VERSIONS_WF(
                            ch_tool_versions,
                            params.publish_dir,
                            params.enable_publish)


    collect_mqc = ch_sentieon_metrics.collect().ifEmpty([])
                        .combine( ch_custom_report.collect().ifEmpty([]))
                        .combine( ch_ado_summary_table.collect().ifEmpty([]))
                        .combine( REPORT_VERSIONS_WF.out.versions.collect().ifEmpty([]))

    params_meta = [
            session_id: workflow.sessionId,
            mode: params.mode,
            genome: params.genome,
            generate_pseudobulk: params.pseudo_bulk_run
    ]

    MULTIQC_WF ( collect_mqc,
                params_meta,
                params.multiqc_config,
                params.publish_dir,
                params.enable_publish
                )

    ch_multiqc_version = MULTIQC_WF.out.version
    MULTIQC_WF.out.report.ifEmpty{ exit 1, "ERROR: cannot generate any MULTIQC Report." }

}

// OnComplete
workflow.onComplete{
    println( "\nPipeline completed successfully.\n\n" )
    output                              = [:]
    output["pipeline_run_name"]         = workflow.runName
    output["pipeline_name"]             = workflow.manifest.name
    output["pipeline_version"]          = workflow.manifest.version
    output["pipeline_session_id"]       = workflow.sessionId
    output["output"]                    = [:]
    output["output"]["bam"]             = [:]


    if ( !params.is_bam ){
        bam_outfile = file("$params.tmp_dir/bam_files.txt")
        bam_outfile_lines = bam_outfile.readLines()
        for ( bam_line : bam_outfile_lines ) {
            def (sample_name, bam_path) = bam_line.split('\t')
            output["output"]["bam"][sample_name] = [:]
            output["output"]["bam"][sample_name]["bam"] = bam_path
        }
    }
    
    def output_json = JsonOutput.toJson(output)
    def output_json_pretty = JsonOutput.prettyPrint(output_json)
    File outputfile = new File("$params.tmp_dir/output.json")
    outputfile.write(output_json_pretty)
    println(output_json_pretty)
}

// OnError
workflow.onError{
    println( "\nPipeline failed.\n\n" )
    output                          = [:]
    output["pipeline_run_name"]     = workflow.runName
    output["pipeline_name"]         = workflow.manifest.name
    output["pipeline_version"]      = workflow.manifest.version
    output["pipeline_session_id"]   = workflow.sessionId
    output["output"]                = [:]

    
    

    def subject = """\
        [nf-somatic_variant_calling-pipeline] FAILED: ${workflow.runName}
        """

    def msg = """\

        Pipeline execution summary 
        --------------------------------
        Script name       : ${workflow.scriptName ?: '-'}
        Script ID         : ${workflow.scriptId ?: '-'}
        Workflow session  : ${workflow.sessionId}
        Workflow repo     : ${workflow.repository ?: '-' }
        Workflow revision : ${workflow.repository ? "$workflow.revision ($workflow.commitId)" : '-'}
        Workflow profile  : ${workflow.profile ?: '-'}
        Workflow cmdline  : ${workflow.commandLine ?: '-'}
        Nextflow version  : ${workflow.nextflow.version}, build ${workflow.nextflow.build} (${workflow.nextflow.timestamp})
        Error Report      : ${workflow.errorReport}
        """
        .stripIndent()
    
    log.info ( msg )
    
    if ( "${params.email_on_fail}" && workflow.exitStatus != 0 ) {
        sendMail(to: "${params.email_on_fail}", subject: subject, body: msg)
    }
}
