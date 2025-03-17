nextflow.enable.dsl=2
params.timestamp = ""

process MERGE_PROCESSED_VCF {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"

    input:
    tuple val(group), path(input_vcfs)
    path(reference)
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val(group), path("merged_group_${group}.vcf.gz*"),emit: merged_vcf
    tuple val(group), path("df_nv_group_${group}.tsv"), emit: df_nv


    script:
    """
    echo -e "Listing files to merge ...";
    
    find . -type f -name "*.vcf.gz" | sort > list_vcf.txt
    cat list_vcf.txt;
    echo -e "Merging VCF files ...";
    bcftools merge --threads ${task.cpus} -m none --file-list list_vcf.txt -g ${reference}/genome.fa | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa -Oz -o merged_group_${group}.vcf.gz
    echo -e "Indexing merged VCF ...";
    bcftools index --threads ${task.cpus} -t merged_group_${group}.vcf.gz
    
    echo -e "Creating table for NV ...";
    
    bcftools query -l merged_group_${group}.vcf.gz  > columns.txt
    cat columns.txt | tr "\\n" "\\t" | sed "s|\\t\$|\\n|" | sed "s|^|\\t|" > df_nv_group_${group}.tsv
    bcftools query --print-header -f '%CHROM\\_%POS\\_%REF\\_%ALT[\\t%AD]\\n' merged_group_${group}.vcf.gz | tail -n+2 >> df_nv_group_${group}.tsv
    """
}