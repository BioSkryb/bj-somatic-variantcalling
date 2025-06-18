nextflow.enable.dsl=2

params.timestamp = ""

process CUSTOM_DATA_PROCESSING {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"


    input:
    tuple val(sample_name), path(dedup_metrics), path(bam), path(bai), path(metrics)
    val(mode)
    val(publish_dir)
    val(enable_publish)


    output:
    path "*_all_sentieonmetrics.txt", emit: metrics
    path("custom_processing_version.yml"), emit: version
    
    script:
    
    // filename = skip_subsampling ? "${sample_name}" : "${sample_name}_${subsample}" 
    // filename = skip_subsampling ? "${sample_name}" : "${sample_name}" 
    
    """
    #! /bin/bash
    set +u
    . /opt/sentieon/cloud_auth.sh no-op

    cat ${sample_name}.cov_sentieonmetrics.sample_interval_summary | tail -n+2 | sed 's|:|\t|g' | sed 's|-|\t|g' | cut -f1-4 > ${sample_name}_chromosome_read_proportions.tsv
    

    parse_metrics_files.py \
       --wgs_metrics_filename ${sample_name}.dedup.wgsmetricsalgo.sentieonmetrics.txt \
       --hs_metrics_filename ${sample_name}.dedup.hsmetricalgo.sentieonmetrics.txt \
       --dedup_metrics_filename ${sample_name}.dedup_sentieonmetrics.txt \
       --summary_metrics_filename ${sample_name}.dedup.alignmentstat_sentieonmetrics.txt \
       --chrom_proportions_filename ${sample_name}_chromosome_read_proportions.tsv \
       --gc_bias_summary_filename ${sample_name}.dedup.gcbias_summary.sentieonmetrics.txt \
       --insert_size_metrics_filename ${sample_name}.dedup.insertsizemetricalgo.sentieonmetrics.txt \
       --sample_name ${sample_name} \
       --mode ${mode}

    mv ${sample_name}_all_metrics.tsv ${sample_name}_all_sentieonmetrics.txt
    
    echo custom_processing: v0.0.1 > custom_processing_version.yml
    """
}

workflow CUSTOM_DATA_PROCESSING_WF {
    take:
        ch_combined_input
        ch_mode
        ch_publish_dir
        ch_enable_publish
        
    main:
        CUSTOM_DATA_PROCESSING ( 
                                ch_combined_input,
                                ch_mode,
                                ch_publish_dir,
                                ch_enable_publish
                             )
                           
    emit:
        metrics = CUSTOM_DATA_PROCESSING.out.metrics
        version = CUSTOM_DATA_PROCESSING.out.version
        
    
}