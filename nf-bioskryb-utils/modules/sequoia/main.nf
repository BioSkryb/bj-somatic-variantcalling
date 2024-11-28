nextflow.enable.dsl=2
params.timestamp = ""

process SEQUOIA {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    

    
    input:
    tuple val(group), path(tabs)
    path(reference)
    val(cutoff_binomial)
    val(cutoff_rho)
    val( publish_dir )
    val( enable_publish )
  
    output:
    path("Sequoia*"), emit: all
    tuple val(group), path("df_filtered_placed_variants_${group}.tsv"), emit:df
    tuple val(group), path("Mat_NR_${group}.tsv"), path("Mat_NV_${group}.tsv"), emit:mat


    script:
    """
    
    
    cat Tab_NR_* | sed 's/NA/0/g'> Mat_NR_${group}.tsv
    
    cat Tab_NV_* | sed 's/NA/0/g' > Mat_NV_${group}.tsv
    
    wc -l Mat_NR_${group}.tsv
    
    wc -l Mat_NV_${group}.tsv

    
    Rscript /usr/local/bin/build_phylogeny.R --genomeFile ${reference}/genome.fa -v Mat_NV_${group}.tsv -r Mat_NR_${group}.tsv --mpboot_path /usr/local/bin/ -n $task.cpus --snv_rho ${cutoff_rho}  --germline_cutoff ${cutoff_binomial}
    
    ls Patient* | while read file; 
    do
    
        name=`echo \${file} | sed 's/Patient/Sequoia_group_${group}_bino${cutoff_binomial}_rho${cutoff_rho}/'`;
    
        mv \${file} \${name};
        
    done
    
    cat *snv_assigned_to_branches.txt  | cut -f 1,2,3,4 | tail -n+2 | sort -u > df_filtered_placed_variants_${group}.tsv
    
    """
    
}