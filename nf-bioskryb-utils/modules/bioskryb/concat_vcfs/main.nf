nextflow.enable.dsl=2
params.timestamp = ""

process CONCAT_VCFS {
    tag "concat_vcfs"
    publishDir "${publish_dir}_${params.timestamp}/tertiary_analyses/variant_annotation/preprocessing", enabled: "$enable_publish"

    
    input:
    path(vcfs)
    path(index)
    val(publish_dir)
    val(enable_publish)
  
    output:
    tuple val("multisample"), path("multisample.vcf.gz*"), emit: multisample_vcf
    
    script:
    """
    
    bcftools concat -a ${vcfs.join(" ")} -O v | bcftools sort -O v -o multisample.vcf
    bgzip multisample.vcf
    tabix -p vcf multisample.vcf.gz
    
    export BCFTOOLS_VER=\$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    echo bcftools merge: \$BCFTOOLS_VER > bcftools_version.yml
    
    """
}