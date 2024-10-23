nextflow.enable.dsl=2
params.timestamp = ""

process SENTIEON_DRIVER_METRICS {
    tag "${sample_name}"
    publishDir "${params.publish_dir}_${params.timestamp}/secondary_analyses/metrics/${sample_name}/alignment_stats", enabled:"$enable_publish"

    input:
    tuple val(sample_name), path(bam), path(bai), path(recal_table_file)
    path fasta_ref
    path base_metrics_intervals
    path wgs_or_target_intervals
    val mode
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val(sample_name), file("*sentieonmetrics*"), emit: metrics_tuple
    path "*sentieonmetrics*", emit: metrics
    tuple val(sample_name), file("*.dedup.alignmentstat_sentieonmetrics.txt"), emit: alignment_metrics
    path("sentiieon_driver_metrics_version.yml"), emit: version
    
    script:
    def bqsr = recal_table_file.name == "dummy_file.txt" ? "" : "-q ${recal_table_file}"
    if (mode == 'exome') {
        """
        set +u
        if [ \$LOCAL != "true" ]; then
            . /opt/sentieon/cloud_auth.sh no-op
        else
            export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER
            echo \$SENTIEON_LICENSE
        fi
    
    
        sentieon driver  -t $task.cpus -r ${fasta_ref}/genome.fa -i ${bam} ${bqsr} \
                --interval ${wgs_or_target_intervals} \
                --algo GCBias --summary ${sample_name}.gcbias_summary.sentieonmetrics.txt ${sample_name}.gcbias.sentieonmetrics.txt \
                --algo AlignmentStat ${sample_name}.dedup.alignmentstat_sentieonmetrics.txt \
                --algo InsertSizeMetricAlgo ${sample_name}.insertsizemetricalgo.sentieonmetrics.txt \
                --algo MeanQualityByCycle ${sample_name}.meanqualitybycycle.sentieonmetrics.txt \
                --algo CoverageMetrics ${sample_name}.cov_sentieonmetrics
                
        #        --algo HsMetricAlgo --targets_list ${wgs_or_target_intervals} --baits_list ${wgs_or_target_intervals} ${sample_name}.hsmetricalgo.sentieonmetrics.txt
        
        # These files are not created in exome. creating dummy files?
        touch ${sample_name}.wgsmetricsalgo.sentieonmetrics.txt
                
                
        export SENTIEON_VER="202308.01"
        echo Sentieon: \$SENTIEON_VER > sentiieon_driver_metrics_version.yml
        """
    } else if (mode == 'wgs') {
        """
        set +u
        if [ \$LOCAL != "true" ]; then
            . /opt/sentieon/cloud_auth.sh no-op
        else
            export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER
            echo \$SENTIEON_LICENSE
        fi
        echo "${bqsr}"
        echo "${recal_table_file.name}"
        
        sentieon driver -t $task.cpus -r ${fasta_ref}/genome.fa -i ${bam} ${bqsr} \
                --interval ${base_metrics_intervals} \
                --algo GCBias --summary ${sample_name}.gcbias_summary.sentieonmetrics.txt ${sample_name}.gcbias.sentieonmetrics.txt \
                --algo AlignmentStat ${sample_name}.dedup.alignmentstat_sentieonmetrics.txt \
                --algo InsertSizeMetricAlgo ${sample_name}.insertsizemetricalgo.sentieonmetrics.txt \
                --algo MeanQualityByCycle ${sample_name}.meanqualitybycycle.sentieonmetrics.txt \
                --algo CoverageMetrics ${sample_name}.cov_sentieonmetrics --omit_base_output --omit_locus_stat --omit_sample_stat
        
        sentieon driver -t $task.cpus -r ${fasta_ref}/genome.fa -i ${bam} ${bqsr} \
                --interval ${wgs_or_target_intervals} \
                --algo WgsMetricsAlgo ${sample_name}.wgsmetricsalgo.sentieonmetrics.txt
                
        # FROM CROMWELL - these files are not created in wgs. creating dummy files?
        # touch ${sample_name}.coveragemetrics.sentieonmetrics.sample_summary
        touch ${sample_name}.hsmetricalgo.sentieonmetrics.txt
        
        export SENTIEON_VER="202308.01"
        echo Sentieon: \$SENTIEON_VER > sentiieon_driver_metrics_version.yml

        """
    } else {
        
        """
        set +u
        
        if [ \$LOCAL != "true" ]; then
            . /opt/sentieon/cloud_auth.sh no-op
        else
            export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER
            echo \$SENTIEON_LICENSE
        fi
        
        sentieon driver  -t $task.cpus -r ${fasta_ref}/genome.fa -i ${bam} \
               --interval ${base_metrics_intervals} \
               --algo WgsMetricsAlgo ${sample_name}.wgsmetricsalgo.sentieonmetrics.txt \
               --algo GCBias --summary ${sample_name}.gcbias_summary.sentieonmetrics.txt ${sample_name}.gcbias.sentieonmetrics.txt \
               --algo AlignmentStat ${sample_name}.dedup.alignmentstat_sentieonmetrics.txt \
               --algo InsertSizeMetricAlgo ${sample_name}.insertsizemetricalgo.sentieonmetrics.txt \
               --algo CoverageMetrics ${sample_name}.cov_sentieonmetrics --omit_base_output --omit_locus_stat --omit_sample_stat
        
        export SENTIEON_VER="202308.01"
        echo Sentieon: \$SENTIEON_VER > sentiieon_driver_metrics_version.yml
        """
    }
}



workflow SENTIEON_DRIVER_METRICS_WF {
    
    take:
        ch_bam
        ch_reference
        ch_base_metrics_intervals
        ch_wgs_or_target_intervals
        ch_mode
        ch_publish_dir
        ch_enable_publish
        
    main:
        SENTIEON_DRIVER_METRICS ( 
                                  ch_bam,
                                  ch_reference,
                                  ch_base_metrics_intervals,
                                  ch_wgs_or_target_intervals,
                                  ch_mode,
                                  ch_publish_dir,
                                  ch_enable_publish
                                )
            
                              
    emit:
        metrics_tuple = SENTIEON_DRIVER_METRICS.out.metrics_tuple
        metrics = SENTIEON_DRIVER_METRICS.out.metrics
        alignment_metrics = SENTIEON_DRIVER_METRICS.out.alignment_metrics
        version = SENTIEON_DRIVER_METRICS.out.version
}

include { CUSTOM_METRICS_MERGE } from '../../../bioskryb/custom_metrics_merge/main.nf' addParams( timestamp: params.timestamp )


workflow {
    
    ch_dummy_file = Channel.fromPath(params.dummy_file, checkIfExists: true).collect()
    
    if (params.bam != "") {
        ch_bam_raw = Channel.fromFilePairs(params.bam, size: -1)
        ch_bam_raw.map{ it -> it.flatten().collect() }
        ch_bam = ch_bam_raw.combine(ch_dummy_file)
    } else if(params.input_csv != "") {
        ch_bam_raw = Channel.fromPath(params.input_csv).splitCsv(header:true)
                                .map { row -> [ row.sampleId, row.bam, row.bam + ".bai" ] }
        ch_bam = ch_bam_raw.combine(ch_dummy_file)
    }
    
    ch_bam.view()
    ch_bam.ifEmpty{ exit 1, "ERROR: No BAM files specified either via --bam or --input_csv" }
    
    SENTIEON_DRIVER_METRICS_WF( 
        ch_bam,
        params.reference,
        params.base_metrics_intervals,
        params.wgs_or_target_intervals,
        params.mode,
        params.publish_dir,
        params.enable_publish
    )
    
    CUSTOM_METRICS_MERGE(
        SENTIEON_DRIVER_METRICS_WF.out.metrics.collect(),
        params.publish_dir,
        params.enable_publish
    )  
}
