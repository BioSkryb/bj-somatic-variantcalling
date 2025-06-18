nextflow.enable.dsl=2
params.timestamp = ""

process SENTIEON_ALIGNMENT {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/alignment", enabled:"$enable_publish"

    input:
    val(genome)
    tuple val(sample_name), path(reads)
    path fasta_ref
    path dbsnp
    path dbsnp_index
    path mills
    path mills_index
    path onekg
    path onekg_index
    val(platform)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(sample_name), path("${sample_name}.bam"), path("${sample_name}.bam.bai"), emit: bam
    tuple val(sample_name), path("*.dedup_sentieonmetrics.txt"), emit: dedup_metrics
    path("*.dedup_sentieonmetrics.txt"), emit: mqc_dedup_metrics
    tuple val(sample_name), path("${sample_name}_recal_data.table"), emit: recal_table
    tuple val(sample_name), path("${sample_name}.bam"), path("${sample_name}.bam.bai"), path("${sample_name}_recal_data.table"), emit: bam_recal_table
    path("sentieon_alignment_version.yml"), emit: version
    
    script:
    if(platform == "Ultima") {
        r1 = "${reads[0]}"
        r2 = "${reads[1]}"
    } else {
        r1 = "${reads[0]}"
        if(reads.size() == 2) {
            r2 = "${reads[1]}"
        } else {
            r2 = ""
        }
    }
    """
    #export PETASUITE_REFPATH=s3://bioskryb-shared-data/dev-resources/petasuite/species
    #export LD_PRELOAD=/usr/lib/petalink.so

    set +u
    if [ \$LOCAL != "true" ]; then
        . /opt/sentieon/cloud_auth.sh no-op
    else
        export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER
        echo \$SENTIEON_LICENSE
    fi
    
    export bwt_max_mem=\$([ ${task.memory.toGiga()} -gt 30 ] && echo "30G" || echo "${task.memory.toGiga()}G")
 
    
    sentieon bwa mem -M -Y -K 2500000000 -R "@RG\\tID:${sample_name}\\tSM:${sample_name}\\tPL:${platform}" -t ${task.cpus} '${fasta_ref}/genome.fa' '${r1}' '${r2}' | sentieon util sort -r '${fasta_ref}/genome.fa' -o '${sample_name}_sorted.bam' -t $task.cpus --sam2bam -i -
    
    sentieon driver -t ${task.cpus} -r '${fasta_ref}/genome.fa' -i ${sample_name}_sorted.bam \
       --algo LocusCollector --fun score_info '${sample_name}.locuscollector_score.gz'

   
    sentieon driver -t ${task.cpus} -r ${fasta_ref}/genome.fa -i ${sample_name}_sorted.bam \
           --algo Dedup --score_info ${sample_name}.locuscollector_score.gz --output_dup_read_name --metrics ${sample_name}.dedup_sentieonmetrics.txt ${sample_name}.tmp_dup_qname.txt.gz
    
    
    sentieon driver -t ${task.cpus} -r ${fasta_ref}/genome.fa -i ${sample_name}_sorted.bam \
           --algo Dedup --rmdup --dup_read_name ${sample_name}.tmp_dup_qname.txt.gz ${sample_name}.bam


    if [[ ${genome} =~ .*GRCh3* ]];
    then
         # run bqsr
        sentieon driver -t ${task.cpus}  -r ${fasta_ref}/genome.fa \
            -i ${sample_name}.bam --algo QualCal \
            -k ${dbsnp} \
            -k ${mills} \
            -k ${onekg} ${sample_name}_recal_data.table
            
         # run bqsr assessment
        sentieon driver -t ${task.cpus}  -r ${fasta_ref}/genome.fa -i ${sample_name}.bam -q ${sample_name}_recal_data.table --algo QualCal \
            -k ${dbsnp} \
            -k ${mills} \
            -k ${onekg} ${sample_name}_recal_data.table.post
    else
                 # run bqsr
        sentieon driver -t ${task.cpus}  -r ${fasta_ref}/genome.fa \
            -i ${sample_name}.bam --algo QualCal \
            ${sample_name}_recal_data.table
            
         # run bqsr assessment
        sentieon driver -t ${task.cpus}  -r ${fasta_ref}/genome.fa -i ${sample_name}.bam -q ${sample_name}_recal_data.table --algo QualCal \
            ${sample_name}_recal_data.table.post
            
    fi
        
    #  plot bqsr outcomes
    sentieon driver -t ${task.cpus} --algo QualCal --plot --before ${sample_name}_recal_data.table --after ${sample_name}_recal_data.table.post ${sample_name}_recal.csv
    sentieon plot QualCal -o ${sample_name}_recal_plots.pdf ${sample_name}_recal.csv

    export SENTIEON_VER="202308.01"
    echo Sentieon: \$SENTIEON_VER > sentieon_alignment_version.yml
    """
}

process SENTIEON_ALGORITHM {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/alignment", enabled:"$enable_publish"

    input:
    val(genome)
    tuple val(sample_name), path(bam), path(bai)
    path fasta_ref
    path dbsnp
    path dbsnp_index
    path mills
    path mills_index
    path onekg
    path onekg_index
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(sample_name), path("${sample_name}_recal_data.table"), emit: recal_table
    tuple val(sample_name), path("${sample_name}.bam"), path("${sample_name}.bam.bai"), path("${sample_name}_recal_data.table"), emit: bam_recal_table
    path("sentieon_alignment_version.yml"), emit: version

    script:
    """
    export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER
    echo \$SENTIEON_LICENSE
    
    if [[ ${genome} =~ .*GRCh3* ]];
    then
         # run bqsr
        sentieon driver -t ${task.cpus}  -r ${fasta_ref}/genome.fa \
            -i ${bam} --algo QualCal \
            -k ${dbsnp} \
            -k ${mills} \
            -k ${onekg} ${sample_name}_recal_data.table
            
         # run bqsr assessment
        sentieon driver -t ${task.cpus}  -r ${fasta_ref}/genome.fa -i ${bam} -q ${sample_name}_recal_data.table --algo QualCal \
            -k ${dbsnp} \
            -k ${mills} \
            -k ${onekg} ${sample_name}_recal_data.table.post
    else
                 # run bqsr
        sentieon driver -t ${task.cpus}  -r ${fasta_ref}/genome.fa \
            -i ${bam} --algo QualCal \
            ${sample_name}_recal_data.table
            
         # run bqsr assessment
        sentieon driver -t ${task.cpus}  -r ${fasta_ref}/genome.fa -i ${bam} -q ${sample_name}_recal_data.table --algo QualCal \
            ${sample_name}_recal_data.table.post
            
    fi
        
    #  plot bqsr outcomes
    sentieon driver -t ${task.cpus} --algo QualCal --plot --before ${sample_name}_recal_data.table --after ${sample_name}_recal_data.table.post ${sample_name}_recal.csv
    sentieon plot QualCal -o ${sample_name}_recal_plots.pdf ${sample_name}_recal.csv
    export SENTIEON_VER="202308.01"
    echo Sentieon: \$SENTIEON_VER > sentieon_alignment_version.yml
    """
}