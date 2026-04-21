nextflow.enable.dsl=2
params.timestamp = ""

process SEQUOIA_SECOND_FILTER {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    input:
    tuple val(group), path(mat_nv), path(mat_nr)
    path(reference)
    val(cutoff_binomial)
    val(cutoff_rho_snp)
    val(cutoff_rho_indel)
    val(min_cov)
    val(max_cov)
    val(gender)
    val(publish_dir)
    val(enable_publish)
  
    output:
    tuple val(group), path("*_filtering_all.txt"), emit: df_filter

    script:
    """

    echo -e "Raw matrices number of lines ...";   
    
    wc -l ${mat_nv}
    
    wc -l ${mat_nr}

    Rscript /usr/local/bin/rscript_4.sequoia_second_pass_filter.R --genomeFile ${reference}/genome.fa -v ${mat_nv} -r ${mat_nr} --mpboot_path /usr/local/bin/ -n $task.cpus --snv_rho ${cutoff_rho_snp} --indel_rho ${cutoff_rho_indel} --germline_cutoff ${cutoff_binomial} --min_cov ${min_cov} --max_cov ${max_cov} --gender ${gender}
    
    ls Patient* | while read file;
    do

        name=`echo \${file} | sed 's/Patient/Sequoia_group_${group}_bino${cutoff_binomial}_rhosnp${cutoff_rho_snp}_rhoindel${cutoff_rho_indel}_mincov${min_cov}_maxcov${max_cov}/'`;

        mv \${file} \${name};

    done

    echo -e "Annotating filtering table with SecondPassFilter column ..."

    filt_file=\$(ls *_filtering_all.txt)
    nr_file=\$(ls *_NR_filtered_all.txt)

    awk 'NR==FNR { if (FNR>1) pass[\$1]=1; next }
         FNR==1  { print \$0, "SecondPassFilter"; next }
         { print \$0, ((\$1 in pass) ? "Pass" : "Fail") }
    ' "\${nr_file}" "\${filt_file}" > filtering_annotated.txt

    mv filtering_annotated.txt "\${filt_file}"

    """
    
}
