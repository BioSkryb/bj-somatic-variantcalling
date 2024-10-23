nextflow.enable.dsl=2
params.timestamp = ""

process PICARD_ADDORREPLACEREADGROUPS {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/", enabled:"$enable_publish"


    input:
    tuple val(sample_name), path(bam), path(bai)
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val(sample_name), path("addorreplacereadgroups/*.bam"), path("addorreplacereadgroups/*.bam.bai"), emit: bam
    path("PICARD_ADDORREPLACEREADGROUPS_version.yml"), emit: version
    
    script:
    """
    #! /bin/bash
    set +u
    
    mkdir addorreplacereadgroups

    picard AddOrReplaceReadGroups -I ${bam} -O addorreplacereadgroups/${sample_name}.bam -RGID ${sample_name} -RGLB lib1 -RGPL Illumina -RGPU unit1 -RGSM ${sample_name}

    cd addorreplacereadgroups
    samtools index -@ $task.cpus ${sample_name}.bam
    
    cd ..

    export SAMTOOLS_VER=\$(samtools --version 2>&1 | sed -e "s/samtools //g")
    echo Samtools: \$SAMTOOLS_VER > PICARD_ADDORREPLACEREADGROUPS_version.yml
    
    export PICARD_VER=\$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    echo Picard: \$PICARD_VER > PICARD_ADDORREPLACEREADGROUPS_version.yml
    """
}

workflow{
    
    if (params.bam != "") {

        ch_bam_raw = Channel.fromFilePairs(params.bam, size: -1)
        ch_bam_raw
                  .map{ it -> it.flatten().collect() }
                  .set{ ch_bam }
    } else if(params.input_csv != "") {
        ch_bam = Channel.fromPath( params.input_csv ).splitCsv( header:true )
                                .map { row -> [ row.sampleId, row.bam, row.bam + ".bai"  ] }
    }
    ch_bam.view()
    ch_bam.ifEmpty{ exit 1, "ERROR: No BAM files specified either via --bam or --input_csv" }

                            
    PICARD_ADDORREPLACEREADGROUPS (
                                    ch_bam,
                                    params.publish_dir,
                                    params.enable_publish
                                )
}
