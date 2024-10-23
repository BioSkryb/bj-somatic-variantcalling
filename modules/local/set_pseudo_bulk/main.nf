params.publish_dir = ""
params.timestamp = ""

process SET_PSEUDO_BULK {
    tag 'SET_PSEUDO_BULK'
    publishDir "${publish_dir}_${params.timestamp}/pseudobulk", enabled:"$enable_publish"

    
    input:
    path(input_csv)
    path(read_counts)
    val(publish_dir)
    val(enable_publish)


    output:
    path "*pseudobulk_input.csv", emit: pseudo_bulk_samples
    path "bulk_bool.txt", emit: bulk_bool

    
    script:
    """
    python3 /scripts/set_pseudobulk_inputs.py -i ${input_csv} -r ${read_counts}
    """
    
}