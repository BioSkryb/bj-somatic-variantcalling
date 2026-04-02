nextflow.enable.dsl=2
params.timestamp = ""

// Per-sample extraction of NR, NV, and GT vectors from an annotated VCF.
// Runs once per sample in parallel; results are grouped downstream and fed to
// CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF which assembles the cross-sample matrices.
//
// Outputs (per sample):
//   <sample>_nr.txt  — one integer per variant: SMPL_PILEUP_NUM_FRAGMENTS_HQ_POSITION
//   <sample>_nv.txt  — one integer per variant: HQ_MQ_BQ_F + HQ_MQ_BQ_R
//   <sample>_gt.txt  — one GT string per variant: [%GT]
//   variant_ids.txt  — CHROM_POS_REF_ALT identifiers (same for every sample in the group;
//                      CREATE_NR_NV_MATRICES uses the copy from the first sample)
process EXTRACT_NR_NV_GT_FROM_ANNOTATED_VCF {
    tag "${group}__${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), val(sample_name), path(vcf), path(tbi)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), val(sample_name),
          path("${sample_name}_nr.txt"),
          path("${sample_name}_nv.txt"),
          path("${sample_name}_gt.txt"),
          path("${sample_name}_variant_ids.txt"),
          emit: per_sample_vectors

    script:
    """
    set -euo pipefail

    # NR: total HQ fragments at position (depth denominator)
    bcftools query \\
        -f '%INFO/SMPL_PILEUP_NUM_FRAGMENTS_HQ_POSITION\\n' \\
        "${vcf}" \\
        | awk '{printf "%d\\n", (\$0+0)}' > "${sample_name}_nr.txt"

    # NV: HQ ALT-supporting fragments (forward + reverse)
    bcftools query \\
        -f '%INFO/SMPL_PILEUP_NUM_FRAGMENTS_HQ_MQ_BQ_F\\t%INFO/SMPL_PILEUP_NUM_FRAGMENTS_HQ_MQ_BQ_R\\n' \\
        "${vcf}" \\
        | awk -F'\\t' '{printf "%d\\n", (\$1+0)+(\$2+0)}' > "${sample_name}_nv.txt"

    # GT: genotype call — used for per-sample presence/absence
    bcftools query \\
        -f '[%GT]\\n' \\
        "${vcf}" > "${sample_name}_gt.txt"

    # Variant IDs — CHROM_POS_REF_ALT (identical across all samples in the group)
    bcftools query \\
        -f '%CHROM\\_%POS\\_%REF\\_%ALT\\n' \\
        "${vcf}" > "${sample_name}_variant_ids.txt"
    """
}
