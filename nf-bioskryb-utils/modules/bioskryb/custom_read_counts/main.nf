params.publish_dir = ""
params.timestamp = ""

process COUNT_READS {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/read_counts", enabled:"$enable_publish"

    input:
    tuple val(sample_name), path(reads)
    val(publish_dir)
    val(enable_publish)

    output:
    path "${sample_name}_read_counts.txt", emit: read_counts

    script:
    """

    echo "Sample names: ${sample_name}"
    echo "Read files: ${reads}"


    python3 /scripts/count_reads.py -s ${sample_name} -r ${reads[0]} -o ${sample_name}_read_counts.txt
    """
}

process COMBINE_READ_COUNTS {
    tag 'COMBINE_READ_COUNTS'
    publishDir "${publish_dir}_${params.timestamp}/read_counts", enabled: "$enable_publish"

    input:
    path(read_counts_files)
    val(publish_dir)
    val(enable_publish)

    output:
    path "combined_read_counts.txt", emit: combined_read_counts

    script:
    """
    # Create the combined file with a single header
    echo "sample_name,reads" > combined_read_counts.txt

    # Append the data from each read_counts file, skipping the header
    for file in ${read_counts_files.flatten().join(' ')}; do
        tail -n +2 \$file >> combined_read_counts.txt
    done
    """
}

process COUNT_READS_FASTQ {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/read_counts", enabled:"$enable_publish"

    input:
    tuple val(sample_name), path(reads)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(sample_name), path(reads), env(total_reads), emit: read_counts
    path "${sample_name}_read_counts.txt", emit: read_counts_file

    script:
    """
    echo "Sample names: ${sample_name}"
    echo "Read files: ${reads}"

    count_reads() {
        if [ -s "\$1" ]; then
            zcat "\$1" 2>/dev/null | awk 'NR%4==1' | wc -l || echo 0
        else
            echo 0
        fi
    }
    
    count_reads "${reads[0]}" > count1.txt &
    pid1=\$!
    count_reads "${reads[1]}" > count2.txt &
    pid2=\$!
    
    wait \$pid1 \$pid2
    
    read1_count=\$(cat count1.txt)
    read2_count=\$(cat count2.txt)
    
    total_reads=\$((read1_count + read2_count))
    if [ "\$total_reads" -eq 0 ]; then
        echo "WARNING: No reads found in ${sample_name}" >&2
    fi

    # Create the output file with header and counts
    echo -e "sample_name,reads" > ${sample_name}_read_counts.txt
    echo -e "${sample_name},\$total_reads" >> ${sample_name}_read_counts.txt
    """
}

workflow CUSTOM_READ_COUNTS_WF{
    take:
        ch_reads
        ch_publish_dir
        ch_enable_publish
        
    main:
        COUNT_READS ( 
                                ch_reads,
                                ch_publish_dir,
                                ch_enable_publish
                             )

        COMBINE_READ_COUNTS ( 
                                COUNT_READS.out.read_counts.collect(),
                                ch_publish_dir,
                                ch_enable_publish
                             )
                           
    emit:
        combined_read_counts = COMBINE_READ_COUNTS.out.combined_read_counts
            
}

workflow COUNT_READS_FASTQ_WF{
    take:
        ch_reads
        ch_publish_dir
        ch_enable_publish
        
    main:
        COUNT_READS_FASTQ ( 
                                ch_reads,
                                ch_publish_dir,
                                ch_enable_publish
                             )

        COMBINE_READ_COUNTS ( 
                                COUNT_READS_FASTQ.out.read_counts_file.collect(),
                                ch_publish_dir,
                                ch_enable_publish
                             )
                           
    emit:
        read_counts = COUNT_READS_FASTQ.out.read_counts
        combined_read_counts = COMBINE_READ_COUNTS.out.combined_read_counts
            
}
