nextflow.enable.dsl=2
params.timestamp = ""

process GATK_VARIANTSTOTABLE {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/tertiary_analyses/variant_annotation/", enabled:"$enable_publish"

    input:
    tuple val(sample_name), path(vcf)
    val(publish_dir)
    val(enable_publish)


    output:
    path("gatk_variantstotable_df.txt"), emit: variants_df
    path("gatk_variantstotable_version.yml"), emit: version
    
    script:
    """

    gatk VariantsToTable -V ${vcf[0]} -F CHROM -F POS -F ID -F REF -F ALT \
                                      -F QUAL -F FILTER -F BaseQRankSum \
                                      -F ClippingRankSum -F ExcessHet \
                                      -F FS -F MQ -F MQRankSum -F QD \
                                      -F ReadPosRankSum -F SOR -F VQSLOD \
                                      -F DP -F AF -F MLEAC -F MLEAF \
                                      -F AN -F AC -F CLNSIG  \
                                      -F CSQ \
                                      -F TYPE -GF AD -GF GT -GF DP -GF PGT -GF PID -GF PL -GF QUAL -O gatk_variantstotable_df.txt

    export GATK_VARIANTSTOTABLE_VER=\$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')      
    echo GATK: \$GATK_VARIANTSTOTABLE_VER > gatk_variantstotable_version.yml
    """
}
