nextflow.enable.dsl=2
params.timestamp = ""

process BULK_GET_VARIANTS_TO_FILTER {
    tag "BULK_GET_VARIANTS_TO_FILTER"
    publishDir "${publish_dir}_${params.timestamp}/BULK_GET_VARIANTS_TO_FILTER/", enabled: "$enable_publish"
    
    input:
    path(input_vcf)
    path(reference)
    val(model_vcf)
    val(publish_dir)
    val(enable_publish)

  
    output:
    path("variants_present_bulk.txt")

    script:
    """

    if [ "${model_vcf}" = "deepvariant" ]; then

        bcftools view --threads ${task.cpus} -i 'GT[*]="alt"' ${input_vcf} | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa | bcftools norm --threads ${task.cpus} -d exact | bcftools view --threads ${task.cpus} -i 'GT[*]="alt"' -Oz -o temp.vcf.gz 
    else

        bcftools view --threads ${task.cpus} -i 'GT[*]="alt"' ${input_vcf} | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa  | bcftools view --threads ${task.cpus} -i 'GT[*]="alt"' -Oz -o temp.vcf.gz

    fi

    bcftools index --threads ${task.cpus} -t temp.vcf.gz

    bcftools query --print-header -f '%CHROM\\_%POS\\_%REF\\_%ALT' temp.vcf.gz | sort -u >  variants_present_bulk.txt


    """
}
