nextflow.enable.dsl=2

params.timestamp = ""

process PUBLISH_INPUT_DATASET {
    tag "publish_input_dataset"
    publishDir "${params.publish_dir}_${params.timestamp}/execution_info", enabled:"$enable_publish"

    input:
        path(file)
        val(publish_dir)
        val(enable_publish)

    output:
        path("input_dataset.csv"), emit: csv_input

    script:
    """
        cat $file > input_dataset.csv
    """
    
}


workflow PUBLISH_INPUT_DATASET_WF{
    take:
        ch_reads_csv
        ch_publish_dir
        ch_enable_publish

    
    main:
        PUBLISH_INPUT_DATASET ( 
                                ch_reads_csv, 
                                ch_publish_dir,
                                ch_enable_publish
                            )
              
    emit:
        csv_input = PUBLISH_INPUT_DATASET.out.csv_input
    
}
