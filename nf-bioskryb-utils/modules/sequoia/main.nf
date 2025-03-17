nextflow.enable.dsl=2
params.timestamp = ""

process SEQUOIA {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"

    input:
    tuple val(group),path(mat_nv), path(mat_nr)
    path(reference)
    val(cutoff_binomial)
    val(cutoff_rho_snp)
    val(cutoff_rho_indel)
    val(min_cov)
    val(max_cov)
    val(publish_dir)
    val(enable_publish)

    output:
    path("Sequoia*"), emit: all
    tuple val(group), path("df_filtered_placed_variants_${group}.tsv"), emit:df
    tuple val(group), path("*both_assigned_to_branches.txt"), path("*both_NV_filtered_all.txt"), path("*both_NR_filtered_all.txt"), path("*_snv_tree_with_branch_length.tree"), emit:bundle_post_vaf_tree

    script:
    """
    
    wc -l ${mat_nv}
    
    wc -l ${mat_nr}
    Rscript /usr/local/bin/sequoia_build_phylogeny.R --genomeFile ${reference}/genome.fa -v ${mat_nv} -r ${mat_nr} --mpboot_path /usr/local/bin/ -n $task.cpus --snv_rho ${cutoff_rho_snp} --indel_rho ${cutoff_rho_indel} --germline_cutoff ${cutoff_binomial} --min_cov ${min_cov} --max_cov ${max_cov}
    
    ls Patient* | while read file; 
    do
    
        name=`echo \${file} | sed 's/Patient/Sequoia_group_${group}_bino${cutoff_binomial}_rhosnp${cutoff_rho_snp}_rhoindel${cutoff_rho_indel}_mincov${min_cov}_maxcov${max_cov}/'`;
    
        mv \${file} \${name};
        
    done
    cat *both_assigned_to_branches.txt  | cut -f 1,2,3,4 | tail -n+2 | sort -u > df_filtered_placed_variants_${group}.tsv
    """

}