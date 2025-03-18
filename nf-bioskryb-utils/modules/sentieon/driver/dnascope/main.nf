nextflow.enable.dsl=2
params.timestamp = ""

process SENTIEON_DNASCOPE {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/variant_calls_dnascope", enabled:"$enable_publish", pattern: "*dnascope.vcf.gz*"


    input:
    tuple val(sample_name), path(bam), path(bai)
    path(fasta_ref)
    path(interval)
    path(dbsnp)
    path(dbsnp_index)
    path(model)
    val(pcrfree)
    val(ploidy)
    val(emit_mode)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(sample_name), path("*_dnascope.vcf.gz*"), emit: vcf
    path("sentieon_dnascope_version.yml"), emit: version


    script:
    dbsnp_param = dbsnp ? '-d ' + dbsnp : ''
    interval_param = interval ? '--interval ' + interval : ''
    //pcrfree = true; change it to pcrfree = false for PCR samples
    pcr_indel_model = pcrfree ? "--pcr_indel_model none" : ''

    """
    set +u
    
    export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER
    echo \$SENTIEON_LICENSE
    
    sentieon driver \
    -t ${task.cpus} \
    -r ${fasta_ref}/genome.fa \
    -i ${bam} \
    ${interval_param} \
    --algo DNAscope \
    ${dbsnp_param} \
    ${pcr_indel_model} \
    --ploidy ${ploidy} \
    --model ${model}/dnascope.model \
    --emit_mode ${emit_mode} \
    TMP_VARIANT.vcf.gz
    
    sentieon driver \
    -t ${task.cpus} \
    -r ${fasta_ref}/genome.fa \
    --algo DNAModelApply \
    --model ${model}/dnascope.model \
    -v TMP_VARIANT.vcf.gz \
    ${sample_name}_dnascope.vcf.gz
        
    export SENTIEON_VER="202308.01"
    echo Sentieon: \$SENTIEON_VER > sentieon_dnascope_version.yml
    """
}


workflow {
    
    log.info """\
      reference                       : ${ params.reference }
      calling_intervals_filename      : ${ params.calling_intervals_filename }
      \n
     """

    if (params.bam != "") {

        ch_bam_raw = Channel.fromFilePairs(params.bam, size: -1)
        ch_bam_raw
                  .map{ it -> it.flatten().collect() }
                  .set{ ch_bam }
    } else if(params.input_csv != "") {
        ch_bam = Channel.fromPath( params.input_csv ).splitCsv( header:true )
                                .map { row -> [ row.biosampleName, row.bam, row.bam + ".bai"  ] }
    }
    ch_bam.view()
    ch_bam.ifEmpty{ exit 1, "ERROR: No BAM files specified either via --bam or --input_csv" }

    SENTIEON_DNASCOPE (    ch_bam,
                     params.reference,
                     params.calling_intervals_filename,
                     params.dbsnp,
                     params.dbsnp_index,
                     params.model,
                     params.pcrfree,
                     params.ploidy,
                     params.dnascope_emit_mode,
                     params.publish_dir,
                     params.enable_publish
                )
}