nextflow.enable.dsl=2
params.timestamp = ""

// Filter a per-sample annotated VCF for mutational signature analysis.
// Retains only variants that pass the pileup quality filters and are not
// homozygous-reference, producing a clean input set for SigProfilerAssignment.
//
// Filter logic:
//   (SMPL_PILEUP_PropClipped_Filter="Pass" || Pileup_Verdict="Pass") && GT != "0/0"
process SUBSET_ANNOTATED_VCFS_FOR_MUTSIG {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", mode: 'copy', enabled: "$enable_publish"

    input:
    tuple val(group), val(sample_name), path(vcf), path(tbi)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(sample_name),
          path("${sample_name}_mutsig.vcf.gz"),
          path("${sample_name}_mutsig.vcf.gz.tbi"),
          emit: filtered_vcf

    script:
    """
    bcftools filter \\
        -i 'INFO/SMPL_PILEUP_PropClipped_Filter="Pass" || INFO/Pileup_Verdict="Pass"' \\
        ${vcf} \\
      | bcftools view \\
        -e 'GT="0/0"' \\
        -O z -o ${sample_name}_mutsig.vcf.gz

    bcftools index -t ${sample_name}_mutsig.vcf.gz
    """
}
