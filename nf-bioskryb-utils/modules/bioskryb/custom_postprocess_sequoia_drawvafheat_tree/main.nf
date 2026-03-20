nextflow.enable.dsl=2
params.timestamp = ""

// Run VAF + digital genotype heatmap postprocessing for SNV, INDEL, and BOTH
// variant types in a single task per group.
// The GT table (sample × variant × genotype) and variant ID list are built once
// from the shared unfiltered NR matrix and reused across all three types.
// Each variant type checks for a valid placement directory independently;
// types with no placement results produce only an empty sentinel directory.
//
// Inputs
//   snv_placement_dirs   — output_snv_placement_* from SEQUOIA_VARIANT_PLACEMENT_SNV
//   indel_placement_dirs — output_indel_placement_* from SEQUOIA_VARIANT_PLACEMENT_INDEL
//   both_placement_dirs  — output_both_placement_* from SEQUOIA_VARIANT_PLACEMENT_BOTH
//   nr_matrix            — NR_annotated_vcf_${group}_unfiltered.tsv (shared)
//   nv_matrix            — NV_annotated_vcf_${group}_unfiltered.tsv (shared)
//   gt_files             — all ${sample}_genotype.tsv for this group (shared)
//
// Outputs (one directory per variant type; always created so output globs resolve)
//   postprocess_snv_${group}/   — vafheatmap.pdf, digitalheatmap.pdf, figures.RDS (or sentinel)
//   postprocess_indel_${group}/ — same
//   postprocess_both_${group}/  — same
process POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group),
          path(snv_placement_dirs),
          path(indel_placement_dirs),
          path(both_placement_dirs),
          path(nr_matrix),
          path(nv_matrix),
          path(gt_files)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group),
          path("postprocess_snv_${group}"),
          path("postprocess_indel_${group}"),
          path("postprocess_both_${group}"),
          emit: postprocess_outputs

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Always create all output directories so output globs always resolve
    for vtype in snv indel both; do
        mkdir -p postprocess_\${vtype}_${group}
        touch postprocess_\${vtype}_${group}/no_results
    done

    # ── Build combined GT table once (shared across all variant types) ────────
    # Source files: \${sample}_genotype.tsv  (header: VARIANT_ID<tab>GT)
    # Extract sample name by stripping the _genotype.tsv suffix.
    for f in *_genotype.tsv; do
        sample="\${f%_genotype.tsv}"
        tail -n +2 "\${f}" | awk -v s="\${sample}" 'BEGIN{OFS="\\t"} {print s, \$1, \$2}'
    done > df_all_gt.tsv

    # Filter GT table to variants present in the unfiltered NR matrix
    awk 'NR>1 {gsub(/"/, "", \$1); print \$1}' ${nr_matrix} > ids.txt
    awk -v OFS="\\t" -v FS="\\t" \\
        'NR==FNR {a[\$0]; next} \$2 in a {print}' \\
        ids.txt df_all_gt.tsv > df_all_gt_chosen.tsv
    echo "GT rows (filtered to NR matrix variants): \$(wc -l < df_all_gt_chosen.tsv)"

    # ── Process each variant type ─────────────────────────────────────────────
    for vtype in snv indel both; do
        echo "=== [\${vtype}] Postprocessing ==="

        # Find the real placement directory (not the sentinel plain file)
        real_dir=\$(ls -d output_\${vtype}_placement_*/ 2>/dev/null | head -1 || true)
        if [ -z "\${real_dir}" ]; then
            echo "[\${vtype}] No placement directory found — skipping"
            continue
        fi

        file_placed=\$(ls "\${real_dir}"/*_assigned_to_branches.txt 2>/dev/null | head -1 || true)
        file_tree=\$(ls "\${real_dir}"/*_tree_with_branch_length.tree 2>/dev/null | head -1 || true)
        if [ -z "\${file_placed}" ] || [ -z "\${file_tree}" ]; then
            echo "[\${vtype}] Missing assigned_to_branches or tree_with_branch_length — skipping"
            continue
        fi

        n_placed=\$(tail -n +2 "\${file_placed}" | wc -l)
        echo "[\${vtype}] Placed variants: \${n_placed}"
        if [ "\${n_placed}" -eq 0 ]; then
            echo "[\${vtype}] No placed variants — skipping"
            continue
        fi

        Rscript /usr/local/bin/rscript_4.postprocess_vaf_drawtree_heatmap_intnodes.R \\
            "\${file_placed}" \\
            "${nv_matrix}" \\
            "${nr_matrix}" \\
            "\${file_tree}" \\
            df_all_gt_chosen.tsv

        mv res_composition.pdf \\
            postprocess_\${vtype}_${group}/res_composition_tree_vafheatmap_${group}_\${vtype}.pdf
        mv res_composition_digital.pdf \\
            postprocess_\${vtype}_${group}/res_composition_tree_digitalheatmap_${group}_\${vtype}.pdf
        mv res_figures.RDS \\
            postprocess_\${vtype}_${group}/res_figures_${group}_\${vtype}.RDS

        echo "[\${vtype}] Done"
    done
    """
}
