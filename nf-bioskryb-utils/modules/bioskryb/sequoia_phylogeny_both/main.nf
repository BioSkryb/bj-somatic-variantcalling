nextflow.enable.dsl=2
params.timestamp = ""

// Runs rscript_sequoia_build_phylogeny_only.R for a single NR/NV matrix pair
// (one filtering-scheme label) using all variants (SNVs + indels).
// Invoked once per label — the subworkflow explodes the matrix channel so each
// label launches an independent compute job for full parallelism.
process SEQUOIA_PHYLOGENY_BOTH {
    tag "${group}_${label}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), val(label), path(nr_matrix), path(nv_matrix)
    val(gender)
    val(vaf_absent)
    val(vaf_present)
    val(create_multi_tree)
    val(mpboot_path)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), path("output_both_${label}"), emit: phylogeny_outputs

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Always create the output dir so the Nextflow output declaration always resolves
    mkdir -p output_both_${label}
    touch output_both_${label}/no_results

    n_vars=\$(tail -n +2 "${nr_matrix}" | wc -l)
    echo "[both:${label}] Total variants: \${n_vars}"

    if [ "\${n_vars}" -eq 0 ]; then
        echo "[both:${label}] No variants — skipping"
        exit 0
    fi

    Rscript /usr/local/bin/rscript_sequoia_build_phylogeny_only.R \\
        --donor_id        "${group}_${label}" \\
        --input_nr        "${nr_matrix}" \\
        --input_nv        "${nv_matrix}" \\
        --output_dir      output_both_${label}/ \\
        --only_snvs       FALSE \\
        --gender          "${gender}" \\
        --vaf_absent      ${vaf_absent} \\
        --vaf_present     ${vaf_present} \\
        --create_multi_tree           ${create_multi_tree} \\
        --mpboot_path     "${mpboot_path}"
    """
}
