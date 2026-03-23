nextflow.enable.dsl=2
params.timestamp = ""

// One task per sample: extract three per-sample VCFs from the group VCF annotated
// by IDENTIFY_GERMLINE_FROM_STATS, each retaining a different germline call-set:
//
//   *_hcgermline_stats.vcf.gz  — HIGH_CONFIDENCE_GERMLINE_FROM_STATS = "Yes"
//                                (Binom stats filter + > 80 % het genotype prevalence)
//
//   *_hcgermline_vep.vcf.gz   — VEP_GERMLINE_HIGH_CONFIDENCE = "Yes"
//                                (VEP AF_filter variants with >= 80 % het prevalence)
//
//   *_hcgermline_bulk.vcf.gz  — RemainingAfterBulk_HIGH_CONFIDENCE = "Yes"
//                                (bulk-derived variants confirmed by >= 80 % het prevalence across cells)
process SUBSET_MERGED_VCF_HIGH_CONFIDENCE_GERMLINE_FROM_STATS {
    tag "${group}:${sample_name}"
    container 'quay.io/biocontainers/bcftools:1.14--h88f3f91_0'
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), val(sample_name), path(annotated_vcf), path(tbi)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), val(sample_name),
          path("${sample_name}_hcgermline_stats.vcf.gz"),
          path("${sample_name}_hcgermline_stats.vcf.gz.tbi"),
          emit: stats_vcf
    tuple val(group), val(sample_name),
          path("${sample_name}_hcgermline_vep.vcf.gz"),
          path("${sample_name}_hcgermline_vep.vcf.gz.tbi"),
          emit: vep_vcf
    tuple val(group), val(sample_name),
          path("${sample_name}_hcgermline_bulk.vcf.gz"),
          path("${sample_name}_hcgermline_bulk.vcf.gz.tbi"),
          emit: bulk_vcf

    script:
    """
    set -euo pipefail

    # ── HIGH_CONFIDENCE_GERMLINE_FROM_STATS = "Yes" ───────────────────────────
    bcftools view -s ${sample_name} ${annotated_vcf} \
      | bcftools filter -i 'INFO/HIGH_CONFIDENCE_GERMLINE_FROM_STATS="Yes"' \
      | bcftools view -Oz -o ${sample_name}_hcgermline_stats.vcf.gz
    bcftools index -t ${sample_name}_hcgermline_stats.vcf.gz

    # ── VEP_GERMLINE_HIGH_CONFIDENCE = "Yes" ─────────────────────────────────
    bcftools view -s ${sample_name} ${annotated_vcf} \
      | bcftools filter -i 'INFO/VEP_GERMLINE_HIGH_CONFIDENCE="Yes"' \
      | bcftools view -Oz -o ${sample_name}_hcgermline_vep.vcf.gz
    bcftools index -t ${sample_name}_hcgermline_vep.vcf.gz

    # ── RemainingAfterBulk_HIGH_CONFIDENCE = "Yes" (bulk-confirmed germline) ──
    bcftools view -s ${sample_name} ${annotated_vcf} \
      | bcftools filter -i 'INFO/RemainingAfterBulk_HIGH_CONFIDENCE="Yes"' \
      | bcftools view -Oz -o ${sample_name}_hcgermline_bulk.vcf.gz
    bcftools index -t ${sample_name}_hcgermline_bulk.vcf.gz

    echo "[SUBSET_MERGED_VCF_HCGERMLINE] ${sample_name} stats : \$(bcftools view -H ${sample_name}_hcgermline_stats.vcf.gz | wc -l) variants"
    echo "[SUBSET_MERGED_VCF_HCGERMLINE] ${sample_name} vep   : \$(bcftools view -H ${sample_name}_hcgermline_vep.vcf.gz   | wc -l) variants"
    echo "[SUBSET_MERGED_VCF_HCGERMLINE] ${sample_name} bulk  : \$(bcftools view -H ${sample_name}_hcgermline_bulk.vcf.gz  | wc -l) variants"
    """
}
