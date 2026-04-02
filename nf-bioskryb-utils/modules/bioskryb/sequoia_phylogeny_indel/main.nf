nextflow.enable.dsl=2
params.timestamp = ""

// Runs rscript_sequoia_build_phylogeny_only.R for a single NR/NV matrix pair
// (one filtering-scheme label) after subsetting to indel rows.
// Invoked once per label — the subworkflow explodes the matrix channel so each
// label launches an independent compute job for full parallelism.
process SEQUOIA_PHYLOGENY_INDEL {
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
    tuple val(group), path("output_indel_${label}"), emit: phylogeny_outputs

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Always create the output dir so the Nextflow output declaration always resolves
    mkdir -p output_indel_${label}
    touch output_indel_${label}/no_results

    # ── Subset to indel rows (len(REF)>1 || len(ALT)>1) ─────────────────────
    awk 'NR==1 { print; next }
         { n=split(\$1,a,"_"); ref=a[n-1]; alt=a[n];
           if (length(ref)>1 || length(alt)>1) print }' "${nr_matrix}" > indel_nr.tsv
    awk 'NR==1 { print; next }
         { n=split(\$1,a,"_"); ref=a[n-1]; alt=a[n];
           if (length(ref)>1 || length(alt)>1) print }' "${nv_matrix}" > indel_nv.tsv

    n_indel=\$(tail -n +2 indel_nr.tsv | wc -l)
    echo "[indel:${label}] Indel variants: \${n_indel}"

    if [ "\${n_indel}" -eq 0 ]; then
        echo "[indel:${label}] No indels — skipping"
        exit 0
    fi

    Rscript /usr/local/bin/rscript_sequoia_build_phylogeny_only.R \\
        --donor_id        "${group}_${label}" \\
        --input_nr        indel_nr.tsv \\
        --input_nv        indel_nv.tsv \\
        --output_dir      output_indel_${label}/ \\
        --only_snvs       FALSE \\
        --gender          "${gender}" \\
        --vaf_absent      ${vaf_absent} \\
        --vaf_present     ${vaf_present} \\
        --create_multi_tree           ${create_multi_tree} \\
        --mpboot_path     "${mpboot_path}"
    """
}
