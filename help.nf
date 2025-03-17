nextflow.enable.dsl=2

def printHeader() {
  
  log.info """\
  BJ-SomaticVariantCalling   P I P E L I N E
  ===================================
  fastq files    : ${ params.reads }
  genome         : ${ params.genome }
  timestamp      : ${ params.timestamp }
  publish_dir    : ${ params.publish_dir }
  \n
  """

}

def helpMessage() {

  def yellow = "\033[0;33m"
  def blue = "\033[0;34m"
  def white = "\033[0;37m"

  log.info """\

    ${blue}BJ-SomaticVariantCalling${white}

    ${yellow}Usage:${white}
        nextflow run main.nf [options]

    ${yellow}Script Options: see nextflow.config${white}

        ${yellow}[required]${white}

        --input_csv                 FILE    Path to input csv file

        
        --publish_dir               DIR     Path of the output directory

        --genome                    STR     Reference genome to use. Available options - GRCh37, GRCh38
                                            DEFAULT: ${params.genome}

        ${yellow}[optional]${white}

        --variant_workflow_type     STR    This parameter sets the appropriate workflow to run based on the provided samples:
                                            - match_normal: if one normal sample is provided for the set of single/tumor samples.
                                            - panel_normal: if more than one normal sample is provided for the set of single/tumor samples.
                                            - pseudobulk: if no normal samples are provided.
                                            - somatic_heurestic_filter: Uses dnascope to make variant calls and then uses a heuristic filter to keep somatic variants 
                                            and filter out germline. It doesn't require bulk samples.
                                            DEFAULT: ${params.variant_workflow_type}

        --panel_of_normal_vcf       FILE    Path to panel of normal vcf file

        --somatic_variant_caller    STR     To select the Variant Caller. The options are tnscope or tnseq.
                                            DEFAULT: ${params.somatic_variant_caller}

        --multiqc_config            DIR     Path to multiqc config
                                            DEFAULT: ${params.multiqc_config}

        --skip_variant_annotation   BOOL    Whether to skip variant annotation
                                            DEFAULT: ${params.skip_annotation}

        --skip_sigprofile           BOOL    Whether to skip Mutational Signature
                                            DEFAULT: ${params.skip_sigprofile}
      
        --help                      BOOL    Display help message
        
    ${yellow}Output:${white}
        - Publish directory with all the results
        - multiqc/ : This section includes output files containing metrics from various tools to create a MultiQC report.
        - secondary_analysis/ : alignment/ Biosample level output containing aligned reads and index file on subsample reads.
                                metrics/ Metrics output from secondary analyses - Alignment, GC bias, Insert Size, Coverage, 
                                and library complexity metrics.
                                The *-pipeline_all_metrics_mqc.txt contains metrics from the All Metrics section of the MultiQC 
                                report found in BaseJumper.
        - variant_calls_<tnscope / tnseq> / : This section includes output files containing variant calls from the TNScope/ TNSeq variant caller.
        

    ${yellow}

    """.stripIndent()


}

workflow{
  printHeader()
  helpMessage()
}
