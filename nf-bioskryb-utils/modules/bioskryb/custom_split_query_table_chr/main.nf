nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_SPLIT_QUERY_TABLE_CHR {
    tag "CUSTOM_SPLIT_QUERY_TABLE_CHR_${sample_name}_${chr}"
    publishDir "${publish_dir}_${params.timestamp}/CUSTOM_SPLIT_QUERY_TABLE_CHR/", enabled: "$enable_publish"

    input:
    tuple val(sample_name), path(query_table), val(group), val(chr)
    val(num_lines)
    val(publish_dir)
    val(enable_publish)


    output:
    path("list_pos_sample_${sample_name}_chr_${chr}_file_*.txt")

    script:
    """
  
    cat ${query_table}  | grep -P "${chr}\\t" | cut -f1,2 > temp.txt
  
    numlines=`wc -l temp.txt | awk '{print \$1}'`;
 
    touch list_pos_sample_${sample_name}_chr_${chr}_file_xaa.txt;
    if [[ \${numlines} != 0 ]];
    then
    	split -l ${num_lines} temp.txt
    
    	ls x* | while read file;
    	do
    
        	mv \${file} list_pos_sample_${sample_name}_chr_${chr}_file_\${file}.txt
    	done
    fi
    """
}