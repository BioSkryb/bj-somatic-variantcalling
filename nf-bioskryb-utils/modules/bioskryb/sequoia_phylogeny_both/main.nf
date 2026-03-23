nextflow.enable.dsl=2
params.timestamp = ""

// For each NR/NV matrix pair produced by CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF,
// use all variants (SNVs + indels) and run rscript_sequoia_build_phylogeny_only.R.
// Outputs land in output_both_<label>/ for each of the 8 filtering schemes.
process SEQUOIA_PHYLOGENY_BOTH {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), path(nr_matrices), path(nv_matrices)
    val(gender)
    val(vaf_absent)
    val(vaf_present)
    val(tree_mut_pval)
    val(keep_ancestral)
    val(split_trees)
    val(genotype_conv_prob)
    val(min_pval_for_true_somatic)
    val(min_variant_reads_shared)
    val(min_vaf_shared)
    val(create_multi_tree)
    val(mpboot_path)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), path("output_both_*"), emit: phylogeny_outputs

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Sentinel so the output glob always resolves even if every label is skipped
    touch output_both_no_results

    for nr_file in NR_annotated_vcf_${group}_*.tsv; do
        label=\$(echo "\$nr_file" | sed 's/NR_annotated_vcf_${group}_//' | sed 's/\\.tsv//')
        nv_file="NV_annotated_vcf_${group}_\${label}.tsv"

        n_vars=\$(tail -n +2 "\$nr_file" | wc -l)
        echo "[both:\${label}] Total variants: \${n_vars}"

        if [ "\$n_vars" -eq 0 ]; then
            echo "[both:\${label}] No variants — skipping"
            continue
        fi

        Rscript /usr/local/bin/rscript_sequoia_build_phylogeny_only.R \\
            --donor_id        "${group}_\${label}" \\
            --input_nr        "\$nr_file" \\
            --input_nv        "\$nv_file" \\
            --output_dir      "output_both_\${label}/" \\
            --only_snvs       FALSE \\
            --gender          "${gender}" \\
            --vaf_absent      ${vaf_absent} \\
            --vaf_present     ${vaf_present} \\
            --tree_mut_pval   ${tree_mut_pval} \\
            --keep_ancestral  ${keep_ancestral} \\
            --split_trees     ${split_trees} \\
            --genotype_conv_prob          ${genotype_conv_prob} \\
            --min_pval_for_true_somatic   ${min_pval_for_true_somatic} \\
            --min_variant_reads_shared    ${min_variant_reads_shared} \\
            --min_vaf_shared              ${min_vaf_shared} \\
            --create_multi_tree           ${create_multi_tree} \\
            --mpboot_path     "${mpboot_path}"
    done
    """
}
