nextflow.enable.dsl=2
params.timestamp = ""

process SENTIEON_DRIVER_COVERAGEMETRICS {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/metrics/${sample_name}", enabled:"$enable_publish"

    input:
    tuple val(sample_name), path(bam), path(bai), path(recal_table_file)
    path fasta_ref
    path target_intervals
    path refseq
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val(sample_name), file("*sentieonmetrics*"), emit: metrics_tuple
    path "*sentieonmetrics*", emit: metrics
    // tuple val(sample_name), file("*.dedup.alignmentstat_sentieonmetrics.txt"), emit: alignment_metrics
    // path("sentiieon_driver_metrics_version.yml"), emit: version
    
    script:
    def bqsr = recal_table_file.name == "dummy_file.txt" ? "" : "-q ${recal_table_file}"
    def refseq_opt = refseq.name == "dummy_file2.txt" ? "" : "--gene_list ${refseq}"
    """
    set +u
    if [ \$LOCAL != "true" ]; then
        . /opt/sentieon/cloud_auth.sh no-op
    else
        export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER
        echo \$SENTIEON_LICENSE
    fi


    sentieon driver -t $task.cpus -r ${fasta_ref}/genome.fa -i ${bam} ${bqsr} \
            --interval ${target_intervals} \
            --algo CoverageMetrics ${sample_name}.coveragemetrics.sentieonmetrics \
            ${refseq_opt} \
            --omit_base_output \
            --cov_thresh 10 --cov_thresh 5


    export SENTIEON_VER="202308.01"
    echo Sentieon: \$SENTIEON_VER > sentiieon_driver_metrics_version.yml
    """

}

workflow SENTIEON_DRIVER_COVERAGEMETRICS_WF{
    take:
        ch_custom_output
        ch_reference
        ch_intervals
        ch_refseq
        ch_publish_dir
        ch_enable_publish
        
    main:
        SENTIEON_DRIVER_COVERAGEMETRICS ( 
                                  ch_custom_output,
                                  ch_reference,
                                  ch_intervals,
                                  ch_refseq,
                                  ch_publish_dir,
                                  ch_enable_publish
                                )
                              
    emit:
        metrics_tuple = SENTIEON_DRIVER_COVERAGEMETRICS.out.metrics_tuple
        metrics = SENTIEON_DRIVER_COVERAGEMETRICS.out.metrics
}

workflow{
    
    ch_dummy_file = Channel.fromPath("$projectDir/../../../../assets/dummy_file.txt", checkIfExists: true).collect()
    file_refseq = file(params.refseq)
    ch_bam_raw = Channel.fromFilePairs(params.bam, size: -1)
                .map{ it -> it.flatten().collect() }
    
    ch_bam = ch_bam_raw.combine( ch_dummy_file )
    ch_bam.view()

    SENTIEON_DRIVER_COVERAGEMETRICS_WF ( 
                                      ch_bam,
                                      params.reference,
                                      params.intervals,
                                      file_refseq,
                                      params.publish_dir,
                                      params.enable_publish
                                  )
}