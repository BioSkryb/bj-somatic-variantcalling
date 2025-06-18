nextflow.enable.dsl=2
params.timestamp = ""

process SUBSET_VCF_VARIANTS {
    tag "${sample_name}_SUBSET_VCF_VARIANTS"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    
    input:
    tuple val(group), val(sample_name), path(vcf_file), path(file_variants)
    path(reference)
    val(publish_dir)
    val(enable_publish)

  
    output:
    tuple val(sample_name), path("*_somatic_filtered_variants.vcf.gz*"), emit: vcf

    script:
    """
    
    bcftools view --threads $task.cpus -H -O v -R ${file_variants} ${vcf_file[0]} > variants.tsv
    
    cat ${file_variants} | awk -v OFS="_" '{print \$1,\$2,\$3,\$4}' > ids.txt
    
    awk -v FS="\\t" -v OFS="\\t" 'NR == FNR { a[\$0]; next }{mid=\$1"_"\$2"_"\$4"_"\$5; if (mid in a){print \$0}}' ids.txt variants.tsv > chosen_variants.tsv

    bcftools view --threads $task.cpus -h ${vcf_file[0]}  | cat - chosen_variants.tsv | bcftools view --threads $task.cpus -Oz -o ${sample_name}_somatic_filtered_variants.vcf.gz
     
    bcftools index -t ${sample_name}_somatic_filtered_variants.vcf.gz
    
    """ 
    
}