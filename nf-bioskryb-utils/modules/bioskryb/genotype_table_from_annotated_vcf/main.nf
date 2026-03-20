nextflow.enable.dsl=2
params.timestamp = ""

// For each per-sample annotated VCF produced by ANNOTATE_SAMPLE_VCF,
// emit a two-column TSV:
//   VARIANT_ID<tab>GT
// where VARIANT_ID = CHROM_POS_REF_ALT and GT is the called genotype (0/0, 0/1, 1/1).
// One Nextflow task per sample, matching the parallelism of ANNOTATE_SAMPLE_VCF.
process GENOTYPE_TABLE_FROM_ANNOTATED_VCF {
    tag "${group}:${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), val(sample_name),
          path(vcf), path(tbi)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), val(sample_name),
          path("${sample_name}_genotype.tsv"),
          emit: genotype_table

    script:
    """
    set -euo pipefail

    printf 'VARIANT_ID\tGT\n' > ${sample_name}_genotype.tsv
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' ${vcf} \
        | awk 'BEGIN{OFS="\t"} {print \$1"_"\$2"_"\$3"_"\$4, \$5}' \
        >> ${sample_name}_genotype.tsv
    """
}
