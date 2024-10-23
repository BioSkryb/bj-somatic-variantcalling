nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_METRICS_MERGE {
    tag "merge_metrics"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/metrics", enabled:"$enable_publish"
    
    
    input:
    path(metrics)
    val(publish_dir)
    val(enable_publish)


    output:
    path "all_sentieonmetrics.txt", emit: metrics

    script:
    
    """
    #! /bin/bash

    merge_metrics_files.py --metrics_dir . --output all_sentieonmetrics.txt

    """
}

workflow{
    
    ch_metrics = Channel.fromPath(params.metrics_dir)
    
    CUSTOM_DATA_PROCESSING_WF(
                                ch_metrics,
                                params.publish_dir,
                                params.enable_publish
                             )   
}

