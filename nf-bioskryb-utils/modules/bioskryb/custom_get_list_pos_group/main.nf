nextflow.enable.dsl=2
params.timestamp = ""

process CUSTOM_GET_LIST_POS_GROUP {
    tag "CUSTOM_GET_LIST_POS_GROUP_${group}_${chr}"
    publishDir "${publish_dir}_${params.timestamp}/CUSTOM_GET_LIST_POS_GROUP/", enabled: "$enable_publish"
    
    input:
    tuple val(group), path(res_tables), val(chr)
    val(num_lines)
    val(publish_dir)
    val(enable_publish)

  
    output:
    path("list_pos_*.txt")

    
    script:
    """
    
    ls 
  
    cat res_table_filtered* |  grep -P "\\t${chr}\\t" | cut -f2,3 | sort -u > temp.txt
    
    numlines=`wc -l temp.txt | awk '{print \$1}'`;

    touch list_pos_group_${group}_chr_${chr}_file_xaa.txt;

    if [[ \${numlines} != 0 ]];
    then

	    split -l ${num_lines} temp.txt
    
    	    ls x* | while read file;
    	    do
    
        	mv \${file} list_pos_group_${group}_chr_${chr}_file_\${file}.txt
        
    	    done
    fi



    """
}
