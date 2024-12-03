nextflow.enable.dsl=2
params.timestamp = ""

process PREPROCESS_VCF {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/PREPROCESS_VCF/", enabled: "$enable_publish"
    
    input:
    tuple val(sample_name), path(input_vcf), val(group)
    path(reference)
    val(publish_dir)
    val(enable_publish)

  
    output:
    tuple val(sample_name), path("${sample_name}_noref_norm_qualinfo.vcf.gz*"), val(group), emit: vcf
    tuple val(sample_name), path("df_query_${sample_name}.tsv"), val(group), emit: query_table

    
    script:
    """
    # filter wildtypes and normalize
    bcftools view --threads ${task.cpus} -i 'GT!="0/0"' ${input_vcf[0]} | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa  | bcftools view --threads ${task.cpus} -i 'GT!="0/0"' | bcftools +fill-tags | bcftools view -Oz -o _temp_${sample_name}.vcf.gz

    bcftools index -t _temp_${sample_name}.vcf.gz

    # add quality information
    #bcftools query -f'%CHROM\\t%POS\\t%QUAL\\n' _temp_${sample_name}.vcf.gz | bgzip -c > ${sample_name}_qual.txt.gz
    #tabix -s1 -b2 -e2 ${sample_name}_qual.txt.gz

    #echo '##FORMAT=<ID=QUAL,Number=1,Type=Float,Description="Per-sample QUAL">' > hdr.txt
    #bcftools annotate --threads ${task.cpus} -a ${sample_name}_qual.txt.gz -c CHROM,POS,FORMAT/QUAL -h hdr.txt -Oz -o ${sample_name}_noref_norm_qualinfo.vcf.gz _temp_${sample_name}.vcf.gz

    #bcftools index -t ${sample_name}_noref_norm_qualinfo.vcf.gz

    mv _temp_${sample_name}.vcf.gz ${sample_name}_noref_norm_qualinfo.vcf.gz
    mv _temp_${sample_name}.vcf.gz.tbi ${sample_name}_noref_norm_qualinfo.vcf.gz.tbi
    
    #Extract only snps and create a table
    bcftools view --types snps  ${sample_name}_noref_norm_qualinfo.vcf.gz  | bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' | tail -n+2 > df_query_${sample_name}.tsv



    """
}
