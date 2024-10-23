nextflow.enable.dsl=2
params.timestamp = ""

process VARIANTQC {
    publishDir "${params.publish_dir}_${params.timestamp}/tertiary_analyses/variant_annotation/variantqc/", enabled:"$enable_publish"

    input:
    tuple val(sample_name), path(vcf)
    path(reference)
    val(publish_dir)
    val(enable_publish)

    output:
    path("output.html"), emit: report
    path("variantqc_version.yml"), emit: version

    script:
    """
    java -jar /DISCVRSeq.jar VariantQC \
     --threads ${task.cpus} \
     -R '${reference}/genome.fa' \
     -V ${vcf[0]} \
     -O output.html
        
    export VARIANTQC_VER=1.20
    echo VariantQC: \$VARIANTQC_VER > variantqc_version.yml
    """
    
}

workflow {
    ch_vcf = Channel.fromFilePairs( params.vcf )
    ch_vcf.view()
    VARIANTQC ( ch_vcf, 
                params.reference, 
                params.publish_dir, 
                params.enable_publish
              )
}
