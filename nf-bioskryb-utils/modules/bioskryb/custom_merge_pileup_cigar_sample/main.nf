nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_MERGE_PILEUP_CIGAR_SAMPLE {
    tag "${sample_name}_${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"

    input:
    tuple val(sample_name), val(group), path(pileup_file), path(cigars_file)
    val(num_lines)
    val(publish_dir)
    val(enable_publish)

    output:
    path("res_grouplevel_pileup_with_cigar*")

    script:
    """
    ls;
    cat pileup_* > df_pileup.tsv
    cat df_name_flag_cigar_* > df_cigars.tsv
    
    awk -v OFS="\\t" 'NR == FNR {  a[\$1]=\$2; next }{n=split(\$8,ids,",");split(\$9,flags,",");outstring=a[ids[1]"_"flags[1]];for(i=2; i<=n; i++){cigar_string=a[ids[i]"_"flags[i]];outstring=outstring","cigar_string;}print \$0,outstring;}' df_cigars.tsv df_pileup.tsv > temp.txt
    
    numlines=`wc -l temp.txt | awk '{print \$1}'`;
    if [[ \${numlines} != 0 ]];
    then
	    split -l ${num_lines} temp.txt
    
    	ls x* | while read file;
    	do
    
        mv \${file} res_grouplevel_pileup_with_cigar_group_${group}_sample_${sample_name}_file_\${file}.tsv
        
    	done
    fi
    """

}