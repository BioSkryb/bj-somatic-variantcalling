nextflow.enable.dsl=2
params.timestamp = ""

// Subset merged VCF to only variants listed in chosen_variants (CHROM_POS_REF_ALT, one per line).
process SUBSET_MERGED_VCF_CHOSEN_VARIANTS {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), path(merged_vcf), path(chosen_variants)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), path("subset_merged_${group}.vcf.gz"), path("subset_merged_${group}.vcf.gz.tbi"), emit: subset_vcf

    script:
    def vcf = merged_vcf[0]
    """
    echo -e "Creating regions from chosen_variants ..."
    awk -F'_' 'NF>=2 {print \$1"\\t"\$2-1"\\t"\$2}' ${chosen_variants} | sort -k1,1 -k2,2n | uniq > regions.bed

    echo -e "Extracting variant lines in regions ..."
    bcftools view --threads ${task.cpus} -R regions.bed -H -O v ${vcf} > variants_body.vcf

    echo -e "Filtering to exact CHROM_POS_REF_ALT match ..."
    awk -v FS="\\t" -v OFS="\\t" 'NR==FNR { a[\$0]; next } { mid=\$1"_"\$2"_"\$4"_"\$5; if (mid in a) print \$0 }' ${chosen_variants} variants_body.vcf > chosen_body.vcf

    echo -e "Building subset VCF ..."
    bcftools view --threads ${task.cpus} -h ${vcf} | cat - chosen_body.vcf | bcftools view --threads ${task.cpus} -Oz -o subset_merged_${group}.vcf.gz

    echo -e "Indexing subset VCF ..."
    bcftools index -t subset_merged_${group}.vcf.gz
    """
}

// Split subset VCF by chromosome for parallel VEP.
process SPLIT_SUBSET_VCF_BY_CHR {
    tag "${group}_${chr}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", mode: 'copy', enabled: "$enable_publish"

    input:
    tuple val(group), path(subset_vcf), path(subset_vcf_tbi), val(chr)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), val(chr), path("subset_${group}_${chr}.vcf.gz"), path("subset_${group}_${chr}.vcf.gz.tbi"), emit: subset_vcf_chr

    script:
    def vcf = subset_vcf[0]
    """
    bcftools view --threads ${task.cpus} -r ${chr} -Oz -o subset_${group}_${chr}.vcf.gz ${vcf}
    bcftools index -t subset_${group}_${chr}.vcf.gz
    """
}

// Run VEP on one chromosome's subset VCF (parallelized by chromosome).
process VEP_ANNOTATE {
    tag "${group}_${chr}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", mode: 'copy', enabled: "$enable_publish"

    input:
    tuple val(group), val(chr), path(vcf_chr), path(vcf_chr_tbi)
    path(reference)
    val(species)
    val(assembly)
    path(cache_dir)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), val(chr), path("vep_${group}_${chr}.vcf.gz"), emit: vep_vcf_chr

    script:
    def vcf_in = vcf_chr[0]
    // cache_dir is staged by Nextflow (e.g. from S3) so it is a local path in the container.
    // It must be the cache root (directory containing the species folder, e.g. .../VEP/).
    def cache_root = cache_dir.toString()
    """
    vep --input_file ${vcf_in} \\
        --output_file vep_${group}_${chr}.vcf.gz \\
        --format vcf \\
        --vcf \\
        --chr ${chr} \\
        --check_existing \\
        --compress_output bgzip \\
        --fork ${task.cpus} \\
        --species ${species} \\
        --assembly ${assembly} \\
        --fasta ${reference}/genome.fa \\
        --cache --dir_cache ${cache_root} \\
        --offline \\
        --af \\
        --max_af \\
        --no_check_variants_order

    # Do not index here; SORT_INDEX_VEP_VCF uses bcftools to sort and index (avoids tbx_index_build3 failures).
    """
}

// Sort and index VEP output with bcftools (avoids tabix tbx_index_build3 failures from transcript-assembly mismatch etc).
process SORT_INDEX_VEP_VCF {
    tag "${group}_${chr}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", mode: 'copy', enabled: "$enable_publish"

    input:
    tuple val(group), val(chr), path(vep_vcf_gz)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), val(chr), path("vep_sorted_${group}_${chr}.vcf.gz"), path("vep_sorted_${group}_${chr}.vcf.gz.tbi"), emit: vep_vcf_chr

    script:
    def vcf_in = vep_vcf_gz[0]
    """
    bcftools sort -Oz -o vep_sorted_${group}_${chr}.vcf.gz ${vcf_in}
    bcftools index -t vep_sorted_${group}_${chr}.vcf.gz
    """
}

// Merge per-chromosome VEP VCFs back into one VCF per group.
process MERGE_VEP_VCF_BY_GROUP {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), path(vep_vcf_list)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), path("group_${group}_vep.vcf.gz"), path("group_${group}_vep.vcf.gz.tbi"), emit: vep_vcf

    script:
    """
    ls vep_sorted_${group}_*.vcf.gz 2>/dev/null | sort -V > vep_list.txt
    bcftools concat -f vep_list.txt -n -Oz -o group_${group}_vep.vcf.gz
    bcftools index -t group_${group}_vep.vcf.gz
    """
}
