nextflow.enable.dsl=2
params.timestamp = ""



process SENTIEON_DRIVER_TNSCOPE_PANEL_OF_NORMAL {
    tag "${bulk_sample_name}_PANEL_OF_NORMAL_TNSCOPE"
    publishDir "${publish_dir}_${params.timestamp}/variant_calls_tnscope/${bulk_sample_name}", enabled:"$enable_publish", pattern: "*.vcf.gz*"

    input:
    tuple val(bulk_sample_name), val(isbulk), path(bulk_bam), path(bulk_recal_table)
    path fasta_ref
    path interval
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val(bulk_sample_name), path("*.vcf.gz*"), emit: vcf
    path("sentieon_tnscope_version.yml"), emit: version
    path("*.vcf.gz"), emit: vcf_only
    path("*.vcf.gz.tbi"), emit: vcf_index_only

    script:

    """
    set +u
    if [ \$LOCAL != "true" ]; then
        . /opt/sentieon/cloud_auth.sh no-op
    else
        export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER
        echo \$SENTIEON_LICENSE
    fi
    
    sentieon driver -t  ${task.cpus} -r ${fasta_ref}/genome.fa \
             -i ${bulk_bam[0]}  -q ${bulk_recal_table} \
             --interval ${interval} \
             --algo TNscope \
             --tumor_sample ${bulk_sample_name} \
             ${bulk_sample_name}_panel_of_normal_tnscope.vcf.gz 
        
    export SENTIEON_VER="202308.01"
    echo sentieon_tnscope: \$SENTIEON_VER > sentieon_tnscope_version.yml
    """
}


process GENERATE_PANEL_OF_NORMAL {
    tag "PANEL_OF_NORMAL"
    publishDir "${publish_dir}_${params.timestamp}/variant_calls_tnscope/panel_of_normal", enabled:"$enable_publish", pattern: "*.vcf.gz*"

    input:
    path(vcf)
    path(vcf_index)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple path("panel_of_normal.vcf.gz"), path("panel_of_normal.vcf.gz.tbi")

    script:
    """
    bcftools merge -m all -f PASS,. --force-samples *.vcf.gz | bcftools plugin fill-AN-AC | bcftools filter -i 'SUM(AC)>1' > panel_of_normal.vcf

    # gzip and index vcf
    bgzip panel_of_normal.vcf
    tabix -p vcf panel_of_normal.vcf.gz
    """
}

process SENTIEON_DRIVER_TNSCOPE_TUMOR_PANEL_OF_NORMAL {
    tag "${sc_sample_name}_PANEL_OF_NORMAL_TNSCOPE"
    publishDir "${publish_dir}_${params.timestamp}/variant_calls_tnscope/${sc_sample_name}_PANEL_OF_NORMAL", enabled:"$enable_publish", pattern: "*.vcf.gz*"

    input:
    tuple val(sc_sample_name), val(isbulk), path(sc_bam), path(sc_recal_table), path(panel_of_normal_vcf), path(panel_of_normal_vcf_index)
    path fasta_ref
    path dbsnp
    path dbsnp_index
    path interval
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val(sc_sample_name), path("*.vcf.gz*"), emit: vcf
    path("sentieon_tnscope_version.yml"), emit: version
    path("*.vcf.gz"), emit: vcf_only
    path("*.vcf.gz.tbi"), emit: vcf_index_only

    script:

    """
    set +u
    if [ \$LOCAL != "true" ]; then
        . /opt/sentieon/cloud_auth.sh no-op
    else
        export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER
        echo \$SENTIEON_LICENSE
    fi
    
    sentieon driver -t  ${task.cpus} -r ${fasta_ref}/genome.fa \
             -i ${sc_bam[0]} -q ${sc_recal_table} \
             --interval ${interval} \
             --algo TNscope \
             --tumor_sample ${sc_sample_name} --pon ${panel_of_normal_vcf} \
             --dbsnp ${dbsnp} ${sc_sample_name}_tnscope.vcf.gz 
    
    rm ${dbsnp}
    rm ${dbsnp_index}
    
    export SENTIEON_VER="202308.01"
    echo sentieon_tnscope: \$SENTIEON_VER > sentieon_tnscope_version.yml
    """
}



process SENTIEON_DRIVER_TNSCOPE {
    tag "${sc_sample_name}_${bulk_sample_name}_TNSCOPE"
    publishDir "${publish_dir}_${params.timestamp}/variant_calls_tnscope/${sc_sample_name}_${bulk_sample_name}", enabled:"$enable_publish", pattern: "*.vcf.gz*"

    input:
    tuple val(group), val(sc_sample_name), path(sc_bam), path(sc_recal_table), val(bulk_sample_name), path(bulk_bam), path(bulk_recal_table)
    path fasta_ref
    path dbsnp
    path dbsnp_index
    path interval
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val(sc_sample_name), path("*.vcf.gz*"), emit: vcf
    path("sentieon_tnscope_version.yml"), emit: version
    path("*.vcf.gz"), emit: vcf_only
    path("*.vcf.gz.tbi"), emit: vcf_index_only

    script:

    """
    set +u
    if [ \$LOCAL != "true" ]; then
        . /opt/sentieon/cloud_auth.sh no-op
    else
        export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER
        echo \$SENTIEON_LICENSE
    fi
    
    sentieon driver -t  ${task.cpus} -r ${fasta_ref}/genome.fa \
             -i ${sc_bam[0]} -q ${sc_recal_table} \
             -i ${bulk_bam[0]}  -q ${bulk_recal_table} \
             --interval ${interval} \
             --algo TNscope \
             --tumor_sample ${sc_sample_name} --normal_sample ${bulk_sample_name} \
             --dbsnp ${dbsnp} ${sc_sample_name}_tnscope.vcf.gz 
    
    rm ${dbsnp}
    rm ${dbsnp_index}
    
    export SENTIEON_VER="202308.01"
    echo sentieon_tnscope: \$SENTIEON_VER > sentieon_tnscope_version.yml
    """
}

// workflow SENTIEON_DRIVER_TNSCOPE_WF{
//     take:
//         ch_combine_outputs
//         ch_reference
//         ch_dbsnp
//         ch_dbsnp_index
//         ch_interval
//         ch_publish_dir
//         ch_enable_publish

//     main:
//         SENTIEON_DRIVER_TNSCOPE ( ch_combine_outputs,
//                                      ch_reference,
//                                      ch_dbsnp,
//                                      ch_dbsnp_index,
//                                      ch_interval,
//                                      ch_publish_dir,
//                                      ch_enable_publish
//                                     )
//     emit:
//         vcf = SENTIEON_DRIVER_TNSCOPE.out.vcf
//         version = SENTIEON_DRIVER_TNSCOPE.out.version
//         vcf_only = SENTIEON_DRIVER_TNSCOPE.out.vcf_only
//         vcf_index_only = SENTIEON_DRIVER_TNSCOPE.out.vcf_index_only

// }


workflow {

    samples_ch = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row -> 
            [ 
                row.biosampleName, 
                row.isbulk,
                [ row.bam, row.bam + ".bai"],
                row.bam.replace(".bam", "_recal_data.table")
            ] 
        }
        .branch { biosampleName, isbulk, bam, recal_table ->
            panel_of_normal: isbulk == "true"
            sc: isbulk == "false"
        }

    if ( params.mode == 'wgs' ) {
            params.calling_intervals_filename   = params.genomes [ params.genome ] [ 'calling_intervals_filename' ]
    }

    if ( params.mode == 'exome' ) {
        params.calling_intervals_filename   = params.genomes [ params.genome ] [ params.exome_panel ] [ 'calling_intervals_filename' ]
    }

    // Create panel of normal if not exist
    if ( params.panel_of_normal_vcf == "") {
        // samples_ch.sc.view()
        samples_ch.panel_of_normal.view()

        SENTIEON_DRIVER_TNSCOPE_PANEL_OF_NORMAL ( samples_ch.panel_of_normal, 
                                                params.reference, 
                                                params.calling_intervals_filename, 
                                                params.publish_dir, 
                                                params.enable_publish ) 

        GENERATE_PANEL_OF_NORMAL ( SENTIEON_DRIVER_TNSCOPE_PANEL_OF_NORMAL.out.vcf_only.collect(), SENTIEON_DRIVER_TNSCOPE_PANEL_OF_NORMAL.out.vcf_index_only.collect(), params.publish_dir, params.enable_publish )

        // combine the tumor to panel of normal pairs
        ch_tumor = samples_ch.sc.combine(GENERATE_PANEL_OF_NORMAL.out)
    
    } else {
        panel_of_normal_ch = Channel.fromPath(params.panel_of_normal_vcf)
                                .map { it -> [ it, it + ".tbi" ] }
        panel_of_normal_ch.view()
        ch_tumor = samples_ch.sc.combine(panel_of_normal_ch)
    }

    ch_tumor.view()

    // run tnscope with generic panel of normal file
    SENTIEON_DRIVER_TNSCOPE_TUMOR_PANEL_OF_NORMAL ( ch_tumor, 
                                 params.reference, 
                                 params.dbsnp, 
                                 params.dbsnp_index, 
                                 params.calling_intervals_filename, 
                                 params.publish_dir, 
                                 params.enable_publish
                                )

}