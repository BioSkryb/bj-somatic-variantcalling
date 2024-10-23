nextflow.enable.dsl=2

params.timestamp = ""

process REPORT_VERSIONS{
    tag "report_versions"
    publishDir "${publish_dir}_${params.timestamp}/execution_info", enabled:"$enable_publish"
    
    
    input:
        path(files)
        val(publish_dir)
        val(enable_publish)
        
    output:
        path("tool_mqc_versions.yml"), emit: versions
        
    script:
    """
        echo ${workflow.manifest.name}: > pipeline_version.yml
        cat pipeline_version.yml ${files} | awk 'NR>1 {print "  " \$0; next} 1' | awk '!seen[\$0]++' > tool_mqc_versions.yml
    """
}

workflow REPORT_VERSIONS_WF{
    take:
        ch_version_files
        ch_publish_dir
        ch_enable_publish
    main:
        REPORT_VERSIONS(
                            ch_version_files,
                            ch_publish_dir,
                            ch_enable_publish
                        )
    emit:
        versions = REPORT_VERSIONS.out.versions
}