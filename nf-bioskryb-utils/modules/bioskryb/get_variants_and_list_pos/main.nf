nextflow.enable.dsl=2
params.timestamp = ""

// Extract all variant IDs (CHROM_POS_REF_ALT, one per line) from merged VCF; multi-allelic expanded.
process GET_VARIANTS_FROM_MERGED_VCF {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), path(merged_vcf)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), path("all_variants_${group}.txt"), emit: all_variants

    script:
    def vcf = merged_vcf instanceof List ? merged_vcf.find { it.name.endsWith('.vcf.gz') } ?: merged_vcf[0] : merged_vcf
    """
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${vcf} | \\
    awk -v OFS="_" 'BEGIN{FS="\\t"} NF>=4 { n=split(\$4,a,","); for(i=1;i<=n;i++) print \$1,\$2,\$3,a[i] }' | \\
    sort -u > all_variants_${group}.txt
    """
}

// From chosen_variants (CHROM_POS_REF_ALT per line), produce list_pos (CHROM\tPOS per line) for pileup.
process GET_LIST_POS_FROM_CHOSEN_VARIANTS {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), path(chosen_variants)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), path("list_pos_variant_${group}.txt"), emit: list_pos

    script:
    """
    awk -F'_' 'NF>=2 {print \$1"\\t"\$2}' ${chosen_variants} | sort -u -k1,1 -k2,2n > list_pos_variant_${group}.txt
    """
}

// Subset df_nv (first column = VariantId) to rows whose VariantId is in chosen_variants.
process FILTER_DF_NV_BY_CHOSEN_VARIANTS {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), path(df_nv), path(chosen_variants)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), path("df_nv_filtered_${group}.tsv"), emit: df_nv_filtered

    script:
    """
    head -n1 ${df_nv} > df_nv_filtered_${group}.tsv
    awk -v FS="\\t" 'NR==FNR { a[\$0]; next } FNR>1 && (\$1 in a) { print \$0 }' ${chosen_variants} ${df_nv} >> df_nv_filtered_${group}.tsv
    """
}
