nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_RSCRIPT_SOMATICSNP_FILTER_3_GROUPLEVEL_PROCESS_PILEUP_SAMPLE {
    tag "${par_file}_${sample_name}_${chr}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    

    
    input:
    tuple val(sample_name), val(chr), val(par_file), path(pileup_file), val(group)
    val( publish_dir )
    val( enable_publish )
  
    output:
    tuple val(group),val(chr), path("res_grouplevel_pileup_chr_${chr}_sample_${sample_name}_file_${par_file}.tsv")

    script:
    """
    

    Rscript /usr/local/bin/rscript_3.grouplevel_process_pileup_sample.R ${pileup_file}
    

    
    mv res.tsv res_grouplevel_pileup_chr_${chr}_sample_${sample_name}_file_${par_file}.tsv

    
    
    """
    
}
