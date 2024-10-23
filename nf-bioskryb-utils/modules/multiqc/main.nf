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
