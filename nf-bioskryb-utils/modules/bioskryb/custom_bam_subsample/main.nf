nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_BAM_SUBSAMPLE {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/alignment/subsample", enabled:"$enable_publish"


    input:
    tuple val(sample_name), path(bam), path(bai), val(subsample)
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val("${sample_name}_${subsample}"), path("*_${subsample}.bam"), path("*_${subsample}.bam.bai"), emit: bam
    tuple val("${sample_name}"), path("*_${subsample}.bam"), path("*_${subsample}.bam.bai"), emit: bam_with_samplename
    path("custom_bam_subsample_version.yml"), emit: version
    
    script:
    
    """
    #! /bin/bash
    set +u
    
    samtools view -f 0x2 --threads $task.cpus ${bam} > samtools_paired.bam
    export PROPERLY_ALIGNED=\$(wc -l samtools_paired.bam | grep -Eo '[0-9]{1,}')
    
    if [ ${subsample} -gt \$PROPERLY_ALIGNED ]
    then
        export PROPORTION=100.9999
    else
	    export PROPORTION=\$(awk -v pf_reads_aligned=\$PROPERLY_ALIGNED -v sample_level=${subsample} 'BEGIN { mratio=sample_level/(pf_reads_aligned); printf(10);printf(mratio,10f) }' )
    fi
    echo "\$PROPORTION | \$PROPERLY_ALIGNED | ${subsample}" 

    samtools view -f 0x2 --threads $task.cpus -s \$PROPORTION ${bam} -b -o ${sample_name}_${subsample}_unsorted.bam  
    samtools sort -l 6 -o ${sample_name}_${subsample}.bam -@ $task.cpus ${sample_name}_${subsample}_unsorted.bam
    samtools index -@ $task.cpus ${sample_name}_${subsample}.bam

    # sentieon util sort -i ${sample_name}_${subsample}_unsorted.bam -o ${sample_name}_${subsample}_sorted.bam
    # sentieon util index ${sample_name}_${subsample}_sorted.bam
    
    export SAMTOOLS_VER=\$(samtools --version 2>&1 |  sed -n -e '1p' | grep -Eo [0-9][.]*[0-9]*)
    echo Samtools: \$SAMTOOLS_VER > custom_bam_subsample_version.yml
    """
}

workflow CUSTOM_BAM_SUBSAMPLE_WF{
    take:
        ch_bam
        ch_publish_dir
        ch_enable_publish
        
    main:
        CUSTOM_BAM_SUBSAMPLE ( 
                                ch_bam,
                                ch_publish_dir,
                                ch_enable_publish
                             )
                           
    emit:
        bam = CUSTOM_BAM_SUBSAMPLE.out.bam
        version = CUSTOM_BAM_SUBSAMPLE.out.version
        
    
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

    
    ch_subsample = Channel.fromList( [params.n_reads] )
    ch_bam_subsample = ch_bam.combine(ch_subsample)
                            
    CUSTOM_BAM_SUBSAMPLE_WF(
                            ch_bam_subsample,
                            params.publish_dir,
                            params.enable_publish
                           )
}
