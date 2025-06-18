nextflow.enable.dsl=2
params.timestamp = ""

process PREPROCESS_VCF {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/PREPROCESS_VCF/", enabled: "$enable_publish"
    
    input:
    tuple val(sample_name), path(input_vcf), val(group)
    path(reference)
    val(model_vcf)
    val(publish_dir)
    val(enable_publish)

  
    output:
    tuple val(sample_name), path("${sample_name}_noref_norm.vcf.gz*"), val(group), emit: vcf
    tuple val(sample_name), path("df_query_${sample_name}.tsv"), val(group), emit: query_table
    tuple val(sample_name), path("df_gt_${sample_name}.tsv"), val(group), emit: df_gt

    
    script:
    """

    if [ "${model_vcf}" = "deepvariant" ]; then

        bcftools view --threads ${task.cpus} -i 'GT[*]="alt"' ${input_vcf[0]} | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa | bcftools norm --threads ${task.cpus} -d exact | bcftools view --threads ${task.cpus} -i 'GT[*]="alt"' -Oz -o temp.vcf.gz 
    else

        bcftools view --threads ${task.cpus} -i 'GT[*]="alt"' ${input_vcf[0]} | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa  | bcftools view --threads ${task.cpus} -i 'GT[*]="alt"' -Oz -o temp.vcf.gz

    fi

    bcftools index --threads ${task.cpus} -t temp.vcf.gz

    echo -e "${sample_name}" > noms.txt;

    bcftools reheader --threads ${task.cpus} -s noms.txt temp.vcf.gz | bcftools view --threads ${task.cpus} -Oz -o ${sample_name}_noref_norm.vcf.gz

    bcftools index --threads ${task.cpus} -t ${sample_name}_noref_norm.vcf.gz

    bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD\\t%DP\\t%GT]\\n' ${sample_name}_noref_norm.vcf.gz | tail -n+2 >> df_nv_nr.tsv

    cat df_nv_nr.tsv  | awk -v OFS="\\t" '{gsub(".*,","",\$5);print \$1"_"\$2"_"\$3"_"\$4,\$5,\$6}'  > df_query_${sample_name}.tsv

    cat df_nv_nr.tsv  | awk -v OFS="\\t" '{print \$1"_"\$2"_"\$3"_"\$4,\$7}' > df_gt_${sample_name}.tsv
    """
}
