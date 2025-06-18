nextflow.enable.dsl=2
params.timestamp = ""

process MULTIQC {
    tag 'multiqc'
    publishDir "${params.publish_dir}_${params.timestamp}/multiqc", enabled:"$enable_publish"

    input:
    path(files)
    val(params_meta)
    path(config)
    val(publish_dir)
    val(enable_publish)


    output:
    path "*multiqc_report.html", emit: report
    path "*_data", emit: data
    path "*_plots", optional:true, emit: plots
    path("multiqc_version.yml"), emit: version
    
    script:

    def params_yamlContent = "report_header_info:\n" + params_meta.collect { key, value -> "  - ${key}: ${value}" }.join("\n")
    """
    cp $config/* .

    # Use the pre-constructed params_yamlContent
    echo -e \"\"\"${params_yamlContent}\"\"\" > user_params.yaml

    multiqc . -n multiqc_report.html -c user_params.yaml

    export MULTIQC_VER=\$(multiqc --version | sed -e "s/multiqc, version //g")
    echo multiqc: \$MULTIQC_VER > multiqc_version.yml
    """
}

process MULTIQC_FOR_PARABRICKS {
    tag 'multiqc'
    publishDir "${params.publish_dir}_${params.timestamp}/multiqc", enabled:"$enable_publish"
    publishDir "${params.report_s3_dir}All_Metrics/wid=${workflow.sessionId}/vc=${params.platform}_${params.mode}/dt=${params.timestamp}/", mode: 'copy', pattern: "all_metrics.parquet", enabled:"$enable_publish"

    cache false

    input:
    path(files)
    val(params_meta)
    path(config)
    val(publish_dir)
    val(enable_publish)


    output:
    path "*multiqc_report.html", emit: report
    path "*_data", emit: data
    path "*_plots", optional:true, emit: plots
    path "all_metrics_mqc.txt", emit: all_metrics
    path "selected_metrics_mqc.txt", emit: selected_metrics
    path("multiqc_version.yml"), emit: version
    path("all_metrics.parquet"), emit: all_metrics_parquet
    
    script:

    def params_yamlContent = "report_header_info:\n" + params_meta.collect { key, value -> "  - ${key}: ${value}" }.join("\n")
    """
    cp $config/* .

    # Use the pre-constructed params_yamlContent
    echo -e \"\"\"${params_yamlContent}\"\"\" > user_params.yaml

    # Call the external Python script
    python3 /scripts/parse_parabricks_metrics.py --base_dir .

    multiqc . -n multiqc_report.html -c user_params.yaml --force

    export MULTIQC_VER=\$(multiqc --version | sed -e "s/multiqc, version //g")
    echo multiqc: \$MULTIQC_VER > multiqc_version.yml
    """
}

workflow MULTIQC_WF{
    take:
        ch_reports
        ch_params_meta
        ch_multiqc_config
        ch_publish_dir
        ch_enable_publish
    main:
        MULTIQC ( 
                  ch_reports,
                  ch_params_meta,
                  ch_multiqc_config,
                  ch_publish_dir,
                  ch_enable_publish
                )
    emit:
        version = MULTIQC.out.version
        data = MULTIQC.out.data
        plots = MULTIQC.out.plots
        report = MULTIQC.out.report
}
