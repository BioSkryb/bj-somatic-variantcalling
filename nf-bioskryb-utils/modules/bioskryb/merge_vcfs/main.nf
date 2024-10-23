nextflow.enable.dsl=2
params.timestamp = ""

process MERGE_VCFS {
    tag "merge_vcfs"
    publishDir "${publish_dir}_${params.timestamp}/tertiary_analyses/variant_annotation/preprocessing", enabled: "$enable_publish"

    
    input:
    path(vcfs)
    path(index)
    val(publish_dir)
    val(enable_publish)
  
    output:
    tuple val("multisample"), path("multisample.vcf.gz*"), emit: multisample_vcf
    
    script:
    int vcf_count = vcfs instanceof List ? vcfs.size() : 1
    """
    # if only one VCF, copy it else merge them    
    if [[ ${vcf_count} -eq 1 ]]; then
        echo "Only one VCF file found, copying it to multisample.vcf"
        gunzip -c ${vcfs[0]} > multisample.vcf
    else
        echo "Merging ${vcf_count} VCF files"
        bcftools merge --threads ${task.cpus} ${vcfs.join(" ")} -O v -o multisample.vcf
    fi

    bgzip multisample.vcf
    tabix -p vcf multisample.vcf.gz
    
    export BCFTOOLS_VER=\$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    echo bcftools merge: \$BCFTOOLS_VER > bcftools_version.yml
    
    """
}