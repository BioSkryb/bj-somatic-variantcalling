nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_BAM_GROUP_PILEUP {
    tag "${sample_name}_${chr}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"

    input:
    tuple val(group), val(sample_name), path(input_bam), path(list_pos), val(chr)
    path(reference)
    val(publish_dir)
    val(enable_publish)


    output:  
    tuple val(group), val(chr), val(sample_name), path("pileup_mq0_bq0_group${group}_${sample_name}_${chr}.txt"), emit:pileup
    tuple val(group), val(chr), val(sample_name), path("${sample_name}_df_nr_group${group}_chr_${chr}.tsv"), emit:df_nr


    script:

    """
    echo -e "Subsetting regions of bam ...";
    date;
    cat ${list_pos} | grep "^${chr}[[:space:]]" | awk -v OFS="\\t" '{print \$1,\$2-1,\$2}' > bed.txt
    samtools view --threads $task.cpus --regions-file bed.txt ${input_bam[0]} -O bam -o regions_subset.bam
    
    samtools index --threads $task.cpus regions_subset.bam
    echo -e "Extracting flags/cigar form subsetted bam ...";
    date;
    samtools view --threads ${task.cpus} regions_subset.bam | cut -f1,2,6 | awk -v OFS="\\t" '{print \$1"_"\$2,\$3}' > df_flags_cigars.tsv
    echo -e "Constructing pileup ...";
    date;
    samtools mpileup regions_subset.bam -f ${reference}/genome.fa --output-BP-5 --output-QNAME --disable-overlap-removal --count-orphans --min-BQ 0 --min-MQ 0 --output-MQ --output-extra FLAG,AS -o pileup_temp.txt -l ${list_pos}
    
    echo -e "Adding cigar information to pileup ...";
    date;
    awk -v OFS="\\t" 'NR == FNR {  a[\$1]=\$2; next }{n=split(\$8,ids,",");split(\$9,flags,",");outstring=a[ids[1]"_"flags[1]];for(i=2; i<=n; i++){cigar_string=a[ids[i]"_"flags[i]];outstring=outstring","cigar_string;}print \$0,outstring;}' df_flags_cigars.tsv pileup_temp.txt > pileup_mq0_bq0_group${group}_${sample_name}_${chr}.txt
    echo -e "Creating df NR ...";
    date;
    cat pileup_mq0_bq0_group${group}_${sample_name}_${chr}.txt | cut -f1,2,4 | awk -v OFS="\\t" '{print \$1"_"\$2,\$3}' > ${sample_name}_df_nr_group${group}_chr_${chr}.tsv
    """
}