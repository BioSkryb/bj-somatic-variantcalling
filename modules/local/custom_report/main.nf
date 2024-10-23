params.publish_dir = ""

params.timestamp = ""

process CUSTOM_REPORT {
    tag 'report'
    label 'process_low'
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/metrics", enabled:"$enable_publish"

    
    input:
    path(files)
    path(input_csv)
    path(ado_full_summary)
    val(publish_dir)
    val(enable_publish)


    output:
    path "*metrics_mqc.txt", emit: mqc
    path("custom_report_version.yml"), emit: version

    
    script:
    """
    python3 /scripts/custom_report.py -s '*sentieonmetrics.txt' -m '*MAPD' -f '*.fastp.json' -i '$input_csv' -a '$ado_full_summary' -p 'wgs' -o 'Summary'

    echo custom_report: v0.7 > custom_report_version.yml
    
    """
    
}

workflow CUSTOM_REPORT_WF{
    take:
        ch_metrics
        ch_input_csv
        ch_ado_full_summary
        ch_publish_dir
        ch_enable_publish
    main:
        CUSTOM_REPORT ( 
                        ch_metrics,
                        ch_input_csv,
                        ch_ado_full_summary,
                        ch_publish_dir,
                        ch_enable_publish
                      )
    emit:
        mqc = CUSTOM_REPORT.out.mqc
        version = CUSTOM_REPORT.out.version
}

workflow{
    
}
