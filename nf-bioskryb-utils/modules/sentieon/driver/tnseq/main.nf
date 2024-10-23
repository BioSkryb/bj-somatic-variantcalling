nextflow.enable.dsl=2
params.timestamp = ""

process SENTIEON_DRIVER_TNSEQ {
    tag "${sc_sample_name}_${bulk_sample_name}_TNSEQ"
    publishDir "${publish_dir}_${params.timestamp}/variant_calls_tnseq/${sc_sample_name}_${bulk_sample_name}", enabled:"$enable_publish", pattern: "*.vcf.gz*"

    input:
    tuple val(group), val(sc_sample_name), path(sc_bam), path(sc_recal_table), val(bulk_sample_name), path(bulk_bam), path(bulk_recal_table)
    path fasta_ref
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val(sc_sample_name), path("*.vcf.gz*"), emit: vcf
    path("sentieon_tnseq_version.yml"), emit: version

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
             --algo TNhaplotyper2 \
             --tumor_sample ${sc_sample_name} --normal_sample ${bulk_sample_name} \
             TMP_OUT_TN_VCF

    sentieon driver -t  ${task.cpus} -r ${fasta_ref}/genome.fa \
            --algo TNfilter \
            --tumor_sample ${sc_sample_name} --normal_sample ${bulk_sample_name} \
            -v TMP_OUT_TN_VCF \
            ${sc_sample_name}_tnseq.vcf.gz
    

    
    export SENTIEON_VER="202308.01"
    echo sentieon_tnseq: \$SENTIEON_VER > sentieon_tnseq_version.yml
    """
}

workflow SENTIEON_DRIVER_TNSEQ_WF{
    take:
        ch_combine_outputs
        ch_reference
        ch_publish_dir
        ch_enable_publish

    main:
        SENTIEON_DRIVER_TNSEQ ( ch_combine_outputs,
                                     ch_reference,
                                     ch_publish_dir,
                                     ch_enable_publish
                                    )
    emit:
        vcf = SENTIEON_DRIVER_TNSEQ.out.vcf
        version = SENTIEON_DRIVER_TNSEQ.out.version

}