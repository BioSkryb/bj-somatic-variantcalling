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
