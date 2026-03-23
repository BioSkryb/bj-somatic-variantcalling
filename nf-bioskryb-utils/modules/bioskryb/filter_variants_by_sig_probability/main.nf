nextflow.enable.dsl=2
params.timestamp = ""

// Filters SNVs in a germline-filtered VCF using artifact-signature probabilities
// from SigProfilerAssignment's Decomposed_MutationType_Probabilities.txt.
//
// Variants whose combined probability mass on the specified artifact signatures
// exceeds max_artifact_prob are removed.  INDELs are always passed through.
//
// Per-active-signature cosine similarities (sample SBS96 spectrum vs COSMIC
// reference vector) are written to a separate long-format TSV for downstream
// cohort-level merging via MERGE_SIG_COSINE_SIMILARITIES.
process FILTER_VARIANTS_BY_SIG_PROBABILITY {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", mode: 'copy', enabled: "$enable_publish"

    input:
    tuple val(sample_name), path(vcf), path(vcf_tbi), path(assignment_dir)
    path(reference)
    val(genome_build)
    val(cosmic_version)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(sample_name), path("${sample_name}_sigfiltered.vcf.gz"),
          path("${sample_name}_sigfiltered.vcf.gz.tbi"),        emit: filtered_vcf
    path("${sample_name}_sig_filter_summary.txt"),              emit: filter_summary
    path("${sample_name}_sig_cosine_similarities.tsv"),         emit: cosine_summary

    script:
    """
    python /scripts/filter_variants_by_signature_probability.py \\
        --vcf                  ${vcf} \\
        --assignment_dir       ${assignment_dir} \\
        --reference            ${reference}/genome.fa \\
        --genome_build         ${genome_build} \\
        --cosmic_version       ${cosmic_version} \\
        --output_vcf           ${sample_name}_sigfiltered.vcf.gz \\
        --output_summary       ${sample_name}_sig_filter_summary.txt \\
        --output_cosine        ${sample_name}_sig_cosine_similarities.tsv
    """
}
