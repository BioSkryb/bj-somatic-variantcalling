nextflow.enable.dsl=2
params.timestamp = ""



process PREPROCESS_VCF {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/tertiary_analyses/mutationalcatalogs/preprocessing", enabled: "$enable_publish"
    
    input:
    tuple val(sample_name), path(input_vcf), path(index)
    path(reference)
    val(publish_dir)
    val(enable_publish)

  
    output:
    tuple val(sample_name), path("${sample_name}_noref.vcf.gz"), path("${sample_name}_noref.vcf.gz.tbi"), emit: vcf

    script:
    """
    # filter wildtypes, include variants with PASS,. and normalize
    if bcftools view -h ${input_vcf[0]} | grep -q '##INFO=<ID=SVTYPE,'; then
        bcftools view --threads ${task.cpus} -i 'GT!="0/0"' -f '.,PASS' ${input_vcf[0]} | bcftools view -e 'SVTYPE="BND" || SVTYPE="INS"' | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa  | bcftools view --threads ${task.cpus} -i 'GT!="0/0"' | bcftools +fill-tags | bcftools view -Oz -o ${sample_name}_noref.vcf.gz
    else
        bcftools view --threads ${task.cpus} -i 'GT!="0/0"' -f '.,PASS' ${input_vcf[0]} | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa  | bcftools view --threads ${task.cpus} -i 'GT!="0/0"' | bcftools +fill-tags | bcftools view -Oz -o ${sample_name}_noref.vcf.gz
    fi
    
    bcftools index -t ${sample_name}_noref.vcf.gz


    """
}

process SIGPROFILERGENERATEMATRIX_MATRIX_GENERATOR {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/tertiary_analyses/mutationalcatalogs", enabled:"$enable_publish"


    input:
    tuple val(sample_name), path(vcf), path(tbi)
    val(genome)
    path(reference)
    path(interval)
    val(publish_dir)
    val(enable_publish)

    output:
    path("*mutationalcatalog.tsv"), emit: mutationalcatalog
    path("sigprofilergeneratematrix_matrix_generator_version.yml"), emit: version


    script:
    interval_param = (interval.toString() != 'dummy_file.txt') ? "--bed_file " + interval : ''
    """
    mkdir vcf
    gunzip -c ${vcf} > vcf/${sample_name}.vcf
    SigProfilerMatrixGenerator matrix_generator mutationalcatalog ${genome} vcf/ ${interval_param} --plot True --volume .
    
    (
    echo -e "Pattern\tMutationType\t${sample_name}"
    for dir in vcf/output/{SBS,DBS,ID}
    do
        if [ -d "\$dir" ]; then
            if ls "\$dir"/*.all 1> /dev/null 2>&1; then
                for file in \$dir/*.all
                do
                    pattern=\$(basename "\$file" .all | sed 's/mutationalcatalog\\.//')
                    awk -v pat="\$pattern" 'BEGIN{FS=OFS="\t"} FNR > 1 {print pat, \$0}' "\$file"
                done
            elif ls "\$dir"/*.region 1> /dev/null 2>&1; then
                for file in \$dir/*.region
                do
                    pattern=\$(basename "\$file" .region | sed 's/mutationalcatalog\\.//')
                    awk -v pat="\$pattern" 'BEGIN{FS=OFS="\t"} FNR > 1 {print pat, \$0}' "\$file"
                done
            else
                echo "No files matching the pattern *.all or *.region found in \$dir."
            fi
        fi
    done
    ) > ${sample_name}_mutationalcatalog.tsv

    export VER=1.2.25
    echo SIGPROFILERGENERATEMATRIX_MATRIX_GENERATOR: \$VER > sigprofilergeneratematrix_matrix_generator_version.yml
    """
}

process CUSTOM_MUTATIONCATALOG_MERGER {
    tag "merge"
    publishDir "${publish_dir}_${params.timestamp}/tertiary_analyses/mutationalcatalogs", enabled:"$enable_publish"
    publishDir "${report_s3_dir}Mutational_Signature/wid=${workflow.sessionId}/vc=${params.variant_caller}_${params.dnascope_model_selection}/dt=${params.timestamp}/", mode: 'copy', pattern: "merged_mutationalcatalog.parquet", enabled:"$enable_publish"

    input:
    path(mutationalcatalogs)
    val(report_s3_dir)
    val(publish_dir)
    val(enable_publish)

    output:
    path("merged_mutationalcatalog.tsv"), emit: merged_mutationalcatalog
    path("merged_mutationalcatalog.parquet"), emit: merged_mutationalcatalog_parquet
    path("merged_mutationalcatalog_mqc.tsv"), emit: merged_mutationalcatalog_mqc
    path("Mutational_Signature_Profile_mqc.png"), emit: sbs96_barplot

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import glob
    import pyarrow as pa
    import pyarrow.parquet as pq

    # Pattern for the files to merge
    file_pattern = '*mutationalcatalog.tsv'

    # Get a list of all files matching the pattern
    files = glob.glob(file_pattern)

    # Initialize an empty DataFrame for the merged data
    merged_df = None

    # Iteratively load each file and merge
    for file in files:
        # Load the file into a DataFrame
        temp_df = pd.read_csv(file, sep='\t')
        
        # Merge with the main DataFrame
        if merged_df is None:
            merged_df = temp_df
        else:
            merged_df = pd.merge(merged_df, temp_df, on=['Pattern', 'MutationType'], how='outer')

    # Fill NaN with 0 if there are missing values
    merged_df.fillna(0, inplace=True)

    # Melt the DataFrame to convert it into the desired format
    df_melted = merged_df.melt(id_vars=['Pattern', 'MutationType'], var_name='SampleName', value_name='Count')

    # Reorder the columns
    df_melted = df_melted[['SampleName', 'Pattern', 'MutationType', 'Count']]

    # Save the converted DataFrame to a new TSV file
    df_melted.to_csv('merged_mutationalcatalog.tsv', sep='\t', index=False)

    df_normalized = df_melted
    
    if df_normalized is None:
        raise Exception("No data was successfully processed")
    
    # Convert Count to float64 (double)
    df_normalized['Count'] = pd.to_numeric(df_normalized['Count'], errors='coerce').astype('float64')
    
    # Calculate totals per sample and mutation pattern
    pattern_sample_totals = df_normalized.groupby(['SampleName', 'Pattern'])['Count'].sum().reset_index()
    pattern_sample_totals.columns = ['SampleName', 'Pattern', 'Pattern_Total']
    
    # Merge totals back to main dataframe
    df_normalized = df_normalized.merge(pattern_sample_totals, on=['SampleName', 'Pattern'])
    
    # Calculate normalized counts per pattern
    df_normalized['Normalized_Count'] = (df_normalized['Count'] / df_normalized['Pattern_Total']) * 1000000
    
    # Calculate percentage within each pattern for each sample
    df_normalized['Percentage'] = df_normalized.groupby(['SampleName', 'Pattern'])['Count'].transform(
        lambda x: (x / x.sum()) * 100
    )
    
    # Round the normalized values
    df_normalized['Normalized_Count'] = df_normalized['Normalized_Count'].round(6)
    df_normalized['Percentage'] = df_normalized['Percentage'].round(6)
    
    # Add metadata columns
    df_normalized['Pattern_Total_Count'] = df_normalized['Pattern_Total']
    
    # Ensure all string columns are string type
    string_columns = ['SampleName', 'Pattern', 'MutationType', 'wid', 'vc', 'dt']
    for col in string_columns:
        if col in df_normalized.columns:
            df_normalized[col] = df_normalized[col].astype(str)
    
    # Ensure numeric columns are float64 (double)
    numeric_columns = ['Count', 'Normalized_Count', 'Percentage', 'Pattern_Total_Count']
    for col in numeric_columns:
        df_normalized[col] = df_normalized[col].astype('float64')
    
    # Create PyArrow schema with explicit types
    fields = [
        ('SampleName', pa.string()),
        ('Pattern', pa.string()),
        ('MutationType', pa.string()),
        ('Count', pa.float64()),
        ('Normalized_Count', pa.float64()),
        ('Percentage', pa.float64()),
        ('Pattern_Total_Count', pa.float64())
    ]
    schema = pa.schema(fields)
    
    # Convert the DataFrame to a PyArrow Table
    table = pa.Table.from_pandas(df_normalized, schema=schema)

    # Save the converted DataFrame to a new Parquet file
    pq.write_table(table, 'merged_mutationalcatalog.parquet')

    # Filter for SBS96 type and save to a text file
    sbs96_df = df_normalized[df_normalized['Pattern'] == 'SBS96']
    sbs96_df = sbs96_df[['MutationType', 'SampleName', 'Pattern', 'Count', 'Pattern_Total', 'Normalized_Count', 'Percentage', 'Pattern_Total_Count']]
    sbs96_df.to_csv('merged_mutationalcatalog_mqc.tsv', sep='\t', index=False)

    # Add sbs96 barplot
    import matplotlib.pyplot as plt

    # Define a color map for the mutation types
    color_map = {
        'C>A': 'skyblue',
        'C>G': 'black',
        'C>T': 'red',
        'T>A': 'grey',
        'T>C': 'lightgreen',
        'T>G': 'lightcoral'
    }

    # Filter the DataFrame for SBS96 pattern
    barplot_data = sbs96_df[sbs96_df['Pattern'] == 'SBS96']

    # Extract the base substitution type from the MutationType
    barplot_data['BaseSubstitution'] = barplot_data['MutationType'].str[2:5]

    # Map the colors to the base substitution types
    barplot_data['Color'] = barplot_data['BaseSubstitution'].map(color_map)

    # Sort the data by the base substitution type within the MutationType
    barplot_data = barplot_data.sort_values(by='BaseSubstitution')

    # Plot the barplot
    plt.figure(figsize=(15, 6))
    plt.bar(barplot_data['MutationType'], barplot_data['Percentage'], color=barplot_data['Color'])
    plt.xlabel('Mutation Type')
    plt.ylabel('Percentage of Single Base Substitutions')
    plt.title('SBS96 Mutation Type Percentages')
    plt.xticks(rotation=90)
    plt.tight_layout()

    # Save the plot
    plt.savefig('Mutational_Signature_Profile_mqc.png')
    """
}

process SIGPROFILE_HEATMAP {

    tag "heatmap"
    publishDir "${publish_dir}_${params.timestamp}/tertiary_analyses/mutationalcatalogs", enabled:"$enable_publish"

    input:
    path(input_data)
    val( publish_dir )
    val( enable_publish )
    val( enough_samples ) // boolean to check if the number of samples is 3 or more

    output:
    path("heatmap_sigprofiler.png")

    when:
    enough_samples

    script:
    """
    
    Rscript /scripts/sigprofile_heatmap.R --data ${input_data}
    """

}

workflow SIGPROFILERGENERATEMATRIX_WF {
    
    take:
        ch_vcf
        ch_interval
        ch_reference
        ch_genome
        ch_sigprofilermatrixgenerator_reference
        ch_report_s3_dir
        ch_publish_dir
        ch_enable_publish
        ch_disable_publish

    main:

        PREPROCESS_VCF (    ch_vcf,
                            ch_reference,
                            ch_publish_dir,
                            ch_disable_publish
                        )

        SIGPROFILERGENERATEMATRIX_MATRIX_GENERATOR (    PREPROCESS_VCF.out.vcf,
                                                        ch_genome,
                                                        ch_sigprofilermatrixgenerator_reference,
                                                        ch_interval,
                                                        ch_publish_dir,
                                                        ch_disable_publish
                                                    )

        CUSTOM_MUTATIONCATALOG_MERGER (    SIGPROFILERGENERATEMATRIX_MATRIX_GENERATOR.out.mutationalcatalog.collect(),
                                                    ch_report_s3_dir,       
                                                    ch_publish_dir,
                                                    ch_enable_publish
                                                )

        // Check if the number of samples is 3 or more
        ch_vcf.count().map { count -> count >= 3 }.set { enough_samples }

        SIGPROFILE_HEATMAP(
            CUSTOM_MUTATIONCATALOG_MERGER.out.merged_mutationalcatalog,
            ch_publish_dir,
            ch_enable_publish,
            enough_samples
        )

        enough_samples.view { enough ->
            if (!enough) {
                log.info "SIGPROFILE_HEATMAP process will not run as there are less than 3 samples."
            }
        }
        
    emit:
        merged_mutationalcatalog = CUSTOM_MUTATIONCATALOG_MERGER.out.merged_mutationalcatalog
        merged_mutationalcatalog_mqc = CUSTOM_MUTATIONCATALOG_MERGER.out.merged_mutationalcatalog_mqc
        sbs96_barplot = CUSTOM_MUTATIONCATALOG_MERGER.out.sbs96_barplot
}


workflow {
    
    log.info """\
      reference                       : ${ params.reference }
      vcf                             : ${ params.vcf }
      \n
     """

    if (params.vcf != "") {

        ch_vcf = Channel.fromFilePairs(params.vcf)
                        .map{ k, v -> [ k, v[0], v[1] ]}
        
    } else if (params.input_csv != "") {
        ch_vcf = Channel.fromPath( params.input_csv ).splitCsv( header:true )
                                .map { row -> [ row.sampleId, row.vcf, row.vcf + ".tbi"  ] }
    }
    ch_vcf.view()
    ch_vcf.ifEmpty{ exit 1, "ERROR: No vcf files specified either via --vcf or --input_csv" }

    ch_interval = Channel.fromPath(params.sigprofilermatrixgenerator_interval, checkIfExists: true).collect()

    SIGPROFILERGENERATEMATRIX_WF (    
            ch_vcf,
            ch_interval,
            params.reference,
            params.genome,
            params.sigprofilermatrixgenerator_reference,
            params.report_s3_dir,
            params.publish_dir,
            params.enable_publish,
            params.disable_publish
    )
}