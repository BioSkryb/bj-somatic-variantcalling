nextflow.enable.dsl=2

// IMPORT MODULES

include { MERGE_VCFS } from '../../modules/bioskryb/merge_vcfs/main.nf' addParams( timestamp: params.timestamp )
include { CONCAT_VCFS } from '../../modules/bioskryb/concat_vcfs/main.nf' addParams( timestamp: params.timestamp )
include { GATK_VARIANTSTOTABLE } from '../../modules/gatk/variantstotable/main.nf' addParams(timestamp: params.timestamp)
// include { CUSTOM_EXPANDED_FIELDS } from '../../modules/bioskryb/custom_expanded_fields/main.nf' addParams(timestamp: params.timestamp)
include { VARIANTQC } from '../../modules/variantqc/main.nf' addParams(timestamp: params.timestamp)

process PREPROCESS_VCF {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/tertiary_analyses/variant_annotation/preprocessing", enabled: "$enable_publish"
    
    input:
    tuple val(sample_name), path(input_vcf)
    path(reference)
    val(publish_dir)
    val(enable_publish)

  
    output:
    tuple val(sample_name), path("${sample_name}_noref_norm_qualinfo.vcf.gz*"), emit: vcf
    path("bcftools_version.yml"), emit: version
    
    script:
    """
    # check if the VCF is for DNAscope or TNscope and preprocess accordingly

    if zgrep -q 'ID=DNAscope' ${input_vcf[0]}; then
        echo "Running DNAscope VCF preprocess"
        # filter wildtypes and normalize
        bcftools view --threads ${task.cpus} -i 'GT!="0/0"' ${input_vcf[0]} | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa | bcftools view --threads ${task.cpus} -i 'GT!="0/0"' | bcftools +fill-tags | bcftools view -Oz -o _temp_${sample_name}.vcf.gz
    elif zgrep -q 'ID=TNscope' ${input_vcf[0]}; then
        echo "Running TNscope VCF preprocess"
        normal_sample_name=\$(bcftools query -l ${input_vcf[0]} | tail -1)
        echo "Normal sample name: \${normal_sample_name}"
        bcftools view --threads ${task.cpus} -s ^\${normal_sample_name} -f PASS -e 'SVTYPE="BND" || SVTYPE="INS"' ${input_vcf[0]} | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa | bcftools +fill-tags | bcftools view -Oz -o _temp_${sample_name}.vcf.gz
    elif zgrep -q 'ID=TNhaplotyper2' ${input_vcf[0]}; then
        echo "Running TNhaplotyper2 VCF preprocess"
        normal_sample_name=\$(bcftools query -l ${input_vcf[0]} | tail -1)
        echo "Normal sample name: \${normal_sample_name}"
        bcftools view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM ${input_vcf[0]} | bcftools view --threads ${task.cpus} -s ^\${normal_sample_name} -f PASS | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa | bcftools +fill-tags | bcftools view -Oz -o _temp_${sample_name}.vcf.gz
    fi

    bcftools index -t _temp_${sample_name}.vcf.gz

    # add quality information
    bcftools query -f'%CHROM\t%POS\t%QUAL\n' _temp_${sample_name}.vcf.gz | bgzip -c > ${sample_name}_qual.txt.gz
    tabix -s1 -b2 -e2 ${sample_name}_qual.txt.gz

    echo '##FORMAT=<ID=QUAL,Number=1,Type=Float,Description="Per-sample QUAL">' > hdr.txt
    bcftools annotate --threads ${task.cpus} -a ${sample_name}_qual.txt.gz -c CHROM,POS,FORMAT/QUAL -h hdr.txt -Oz -o ${sample_name}_noref_norm_qualinfo.vcf.gz _temp_${sample_name}.vcf.gz

    bcftools index -t ${sample_name}_noref_norm_qualinfo.vcf.gz

    export BCFTOOLS_VER=\$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    echo  BCFtools: \$BCFTOOLS_VER > bcftools_version.yml


    """
}

process SPLIT_VCF_BY_CHR {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/tertiary_analyses/variant_annotation/split_vcfs", enabled: "$enable_publish"

    input:
    tuple val(sample_name), path(vcf)
    val(publish_dir)
    val(enable_publish)

    output:
    path("${sample_name}_chr*.vcf.gz"), emit: vcf

    script:
    """
    #!/bin/bash
    set -xe
    chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)

    #bgzip ${vcf}
    #tabix -p vcf ${vcf}.gz

    for chr in "\${chromosomes[@]}"
    do
        bcftools view ${vcf} --regions "\$chr" -Oz -o ${sample_name}_"\$chr".vcf.gz
    done
    """
}

process VARIANT_ANNOTATION_VEP {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/tertiary_analyses/variant_annotation/vep", enabled:"$enable_publish", pattern: "*summary.html"
    
    input:
    tuple val(sample_name), path(vcf)
    path(vep_cache)
    val(publish_dir)
    val(enable_publish)
    
    output:
    path("*variants_annotated.vcf.gz"), emit: vcf
    path("*variants_annotated.vcf.gz.tbi"), emit: index
    path("*_summary.html")              , emit: report
    path("vep_version.yml"), emit: version
    
    script:
    def memory = "${task.memory.giga}g"
    """
    #!/bin/bash
    set -xe

    tabix -p vcf ${vcf}

    vep -i ${vcf} -o ${sample_name}_variants_annotated.vcf.gz --fork ${task.cpus} --cache --dir_cache ${vep_cache} --vcf --compress_output bgzip --format vcf --assembly GRCh38 --force_overwrite --offline --everything

    tabix -p vcf ${sample_name}_variants_annotated.vcf.gz



    export VEP_VER=111.0
    echo VEP: \$VEP_VER > vep_version.yml
    """
}


workflow VARIANT_ANNOTATION_WF {
    take:
        ch_vcf
        ch_genome
        ch_reference
        ch_vep_cache
        ch_publish_dir
        ch_enable_publish
        ch_disable_publish
        

    main:
    
    
        /*
        ========================================================================================
            FILTER WILDTYPES, NORMALIZE, AND ADD QUALITY INFORMATION
        ========================================================================================
        */

        

        PREPROCESS_VCF (
                            ch_vcf,
                            ch_reference,
                            params.publish_dir,
                            params.disable_publish
        )

        ch_vcf = PREPROCESS_VCF.out.vcf
    
        /*
        ========================================================================================
            MERGING AND QC
        ========================================================================================
        */

        vcf_file_channel = Channel.empty()
        vcf_index_channel = Channel.empty()
        results = ch_vcf.multiMap { name, pair -> 
                                                files: pair[0]
                                                index: pair[1]            
                                            }
        MERGE_VCFS (
                    vcf_file_channel.mix(results.files).collect(),
                    vcf_index_channel.mix(results.index).collect(),
                    ch_publish_dir,
                    ch_disable_publish
                )

        
        VARIANTQC (
                MERGE_VCFS.out.multisample_vcf,
                ch_reference,
                ch_publish_dir,
                ch_enable_publish
            )

         /*
        ========================================================================================
            ANNOTATING
        ========================================================================================
        */
        
        if ( ch_genome == 'GRCh38' ) {

            SPLIT_VCF_BY_CHR (
                                MERGE_VCFS.out.multisample_vcf,
                                ch_publish_dir,
                                ch_disable_publish
                            )

            // SPLIT_VCF_BY_CHR.out.vcf.view()

            ch_split_vcfs = SPLIT_VCF_BY_CHR.out.vcf.flatMap{ files -> 
                files.collect{ file ->
                    def chrPart = file.baseName - /^multisample_/
                    tuple(chrPart, file) 
                }
            }
            
            // ch_split_vcfs.view()

            VARIANT_ANNOTATION_VEP (
                        ch_split_vcfs,
                        ch_vep_cache,
                        ch_publish_dir,
                        ch_enable_publish
                     )
            
            CONCAT_VCFS (
                VARIANT_ANNOTATION_VEP.out.vcf.collect(),
                VARIANT_ANNOTATION_VEP.out.index.collect(),
                ch_publish_dir,
                ch_enable_publish
            )

            ch_multisample_vcf = CONCAT_VCFS.out.multisample_vcf
        }

        GATK_VARIANTSTOTABLE ( 
                            ch_multisample_vcf,
                            ch_publish_dir,
                            ch_enable_publish
                         )

        // Masking out the custom expand fields for now as it does not work with the large files.
        // CUSTOM_EXPANDED_FIELDS ( 
        //                     GATK_VARIANTSTOTABLE.out.variants_df,
        //                     ch_publish_dir,
        //                     ch_enable_publish
        //                  )


        // ch_variant_summary_df = CUSTOM_EXPANDED_FIELDS.out.variants_df
        
    emit:
        vcf = ch_multisample_vcf
        // variants_df = ch_variant_summary_df
        report = VARIANT_ANNOTATION_VEP.out.report
        bcftools_version = PREPROCESS_VCF.out.version
        vep_version = VARIANT_ANNOTATION_VEP.out.version
        gatk_version = GATK_VARIANTSTOTABLE.out.version
        variantqc_version = VARIANTQC.out.version
        
}
