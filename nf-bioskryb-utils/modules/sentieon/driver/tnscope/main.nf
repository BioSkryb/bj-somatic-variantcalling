nextflow.enable.dsl=2
params.timestamp = ""

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

workflow SENTIEON_DRIVER_TNSCOPE_WF{
    take:
        ch_combine_outputs
        ch_reference
        ch_dbsnp
        ch_dbsnp_index
        ch_interval
        ch_publish_dir
        ch_enable_publish

    main:
        SENTIEON_DRIVER_TNSCOPE ( ch_combine_outputs,
                                     ch_reference,
                                     ch_dbsnp,
                                     ch_dbsnp_index,
                                     ch_interval,
                                     ch_publish_dir,
                                     ch_enable_publish
                                    )
    emit:
        vcf = SENTIEON_DRIVER_TNSCOPE.out.vcf
        version = SENTIEON_DRIVER_TNSCOPE.out.version
        vcf_only = SENTIEON_DRIVER_TNSCOPE.out.vcf_only
        vcf_index_only = SENTIEON_DRIVER_TNSCOPE.out.vcf_index_only

}