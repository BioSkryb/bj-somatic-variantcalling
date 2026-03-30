nextflow.enable.dsl=2
params.timestamp = ""

// Filter a per-sample annotated VCF for mutational signature analysis.
// Retains only variants that are not homozygous-reference (genotype filter,
// always applied). Optionally also applies a pileup quality filter.
//
// Filter logic (quality_filter = true):
//   Pileup_Verdict="Pass" && GT != "0/0"
// Filter logic (quality_filter = false):
//   GT != "0/0"
process SUBSET_ANNOTATED_VCFS_FOR_MUTSIG {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", mode: 'copy', enabled: "$enable_publish"

    input:
    tuple val(group), val(sample_name), path(vcf), path(tbi)
    val(quality_filter)   // true → apply pileup quality filter; false → genotype filter only
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(sample_name),
          path("${sample_name}_mutsig.vcf.gz"),
          path("${sample_name}_mutsig.vcf.gz.tbi"),
          emit: filtered_vcf

    script:
    if ( quality_filter.toString() == "true" )
    """
    bcftools filter \\
        -i 'INFO/Pileup_Verdict="Pass"' \\
        ${vcf} \\
      | bcftools view \\
        -e 'GT="0/0"' \\
        -O z -o ${sample_name}_mutsig.vcf.gz

    bcftools index -t ${sample_name}_mutsig.vcf.gz
    """
    else
    """
    bcftools view \\
        -e 'GT="0/0"' \\
        ${vcf} \\
        -O z -o ${sample_name}_mutsig.vcf.gz

    bcftools index -t ${sample_name}_mutsig.vcf.gz
    """
}
