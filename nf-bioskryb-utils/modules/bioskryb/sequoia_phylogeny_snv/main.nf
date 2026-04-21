nextflow.enable.dsl=2
params.timestamp = ""

// Runs rscript_sequoia_build_phylogeny_only.R for a single NR/NV matrix pair
// (one filtering-scheme label) after subsetting to SNV rows.
// Invoked once per label — the subworkflow explodes the matrix channel so each
// label launches an independent compute job for full parallelism.
process SEQUOIA_PHYLOGENY_SNV {
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
    tuple val(group), path("output_snv_${label}"), emit: phylogeny_outputs

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Always create the output dir so the Nextflow output declaration always resolves
    mkdir -p output_snv_${label}
    touch output_snv_${label}/no_results

    # ── Subset to SNV rows (len(REF)==1 && len(ALT)==1) ──────────────────────
    awk 'NR==1 { print; next }
         { n=split(\$1,a,"_"); ref=a[n-1]; alt=a[n];
           if (length(ref)==1 && length(alt)==1) print }' "${nr_matrix}" > snv_nr.tsv
    awk 'NR==1 { print; next }
         { n=split(\$1,a,"_"); ref=a[n-1]; alt=a[n];
           if (length(ref)==1 && length(alt)==1) print }' "${nv_matrix}" > snv_nv.tsv

    n_snv=\$(tail -n +2 snv_nr.tsv | wc -l)
    echo "[snv:${label}] SNV variants: \${n_snv}"

    if [ "\${n_snv}" -eq 0 ]; then
        echo "[snv:${label}] No SNVs — skipping"
        exit 0
    fi

    Rscript /usr/local/bin/rscript_sequoia_build_phylogeny_only.R \\
        --donor_id        "${group}_${label}" \\
        --input_nr        snv_nr.tsv \\
        --input_nv        snv_nv.tsv \\
        --output_dir      output_snv_${label}/ \\
        --only_snvs       TRUE \\
        --gender          "${gender}" \\
        --vaf_absent      ${vaf_absent} \\
        --vaf_present     ${vaf_present} \\
        --create_multi_tree           ${create_multi_tree} \\
        --mpboot_path     "${mpboot_path}"
    """
}
