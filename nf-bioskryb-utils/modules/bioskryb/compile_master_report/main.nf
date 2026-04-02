nextflow.enable.dsl=2
params.timestamp = ""

// Compile a single master PDF report for one group by assembling outputs from:
//   1. POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE  — VAF heatmaps (SNV / INDEL / BOTH)
//   2. POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE  — digital genotype heatmaps (SNV / INDEL / BOTH)
//   3. PLOT_SIGNATURE_BARGRAPHS             — combined signature activity bar chart PNG (optional)
//   4. TREES_COMPARE_SIMILARITIES           — tree topology comparison master PDF
//   5. CUSTOM_VARIANT_FILTER_PROVENANCE     — combined variant filter provenance report
//
// Each section is preceded by a navy title page generated with R (base graphics).
// PNGs are converted to PDF via R (png::readPNG + rasterImage) before assembly.
// Final output: master_report_${group}.pdf assembled with pdfunite.
//
// Inputs:
//   postprocess_dirs     — list of postprocess_snv/indel/both_${group}/ directories
//   tree_comparison_pdf  — ${group}_tree_comparison_master.pdf
//   combined_report_pdf  — Patient_filter_report_combined_${group}.pdf
//   zero_filtered_png    — zero_filtered_signature_activities.png  (or /dev/null sentinel)
//   cosine_filtered_png  — cosine_filtered_signature_activities.png (or /dev/null sentinel)
//   prevalence_plot_png  — germline_prevalence_distributions_${group}.png
//   ado_summary_png      — ADO_germline_comparison.png
//   matrix_scheme_pdf    — matrix_scheme_summary_${group}.pdf (or /dev/null sentinel)
process COMPILE_MASTER_REPORT {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group),
          path(postprocess_dirs),
          path(tree_comparison_pdf, stageAs: 'tree_comparison.pdf'),
          path(combined_report_pdf, stageAs: 'combined_report.pdf'),
          path(prevalence_plot_png, stageAs: 'prevalence_plot.png'),
          path(matrix_scheme_pdf,   stageAs: 'matrix_scheme.pdf')
    path(zero_filtered_png,   stageAs: 'zero_filtered_sig.png')
    path(cosine_filtered_png, stageAs: 'cosine_filtered_sig.png')
    path(ado_summary_png,     stageAs: 'ado_summary.png')
    path(variant_filter_plot, stageAs: 'variant_filter_plot.pdf')     // variant_filter_plot_${group}.pdf — placed on page 2
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), path("master_report_${group}.pdf"), emit: master_report

    script:
    """
    set -euo pipefail

    parts=""

    # add_section <title> <output_pdf>
    add_section() {
        local title=\$1 out=\$2
        Rscript /usr/local/bin/make_title.R "\$out" "\$title" "${group}"
        parts="\$parts \$out"
    }

    # add_pdf <path>  — appends only if the file exists and is non-empty
    add_pdf() {
        local f=\$1
        if [ -f "\$f" ] && [ -s "\$f" ]; then
            parts="\$parts \$f"
        else
            echo "[compile_master_report] SKIP (missing/empty): \$f"
        fi
    }

    # ── PAGE 1: Master title ──────────────────────────────────────────────────
    Rscript /usr/local/bin/make_title.R "master_title.pdf" "Somatic SNP/INDEL Filter Analysis Report" "${group}"
    parts="master_title.pdf"

    # ── PAGE 2: Variant Filtering Funnel ─────────────────────────────────────
    add_pdf "${variant_filter_plot}"

    # ── SECTION 0: Matrix Filtering Scheme Summary ───────────────────────────
    add_section "Matrix Filtering Scheme Summary" "sec00_matrix_scheme_summary.pdf"
    add_pdf "${matrix_scheme_pdf}"

    # ── SECTION 1: VAF Heatmaps ───────────────────────────────────────────────
    add_section "VAF Heatmaps - Phylogenetic Trees" "sec01_vaf_heatmaps.pdf"
    for vtype in snv indel both; do
        add_pdf "postprocess_\${vtype}_${group}/res_composition_tree_vafheatmap_${group}_\${vtype}.pdf"
    done

    # ── SECTION 2: Digital Genotype Heatmaps ─────────────────────────────────
    add_section "Digital Genotype Heatmaps - Phylogenetic Trees" "sec02_digital_heatmaps.pdf"
    for vtype in snv indel both; do
        add_pdf "postprocess_\${vtype}_${group}/res_composition_tree_digitalheatmap_${group}_\${vtype}.pdf"
    done

    # ── SECTION 3: Mutational Signature Activities ────────────────────────────
    add_section "Mutational Signature Activities" "sec03_signatures.pdf"
    if [ -s "${zero_filtered_png}" ]; then
        Rscript /usr/local/bin/png_to_pdf.R "${zero_filtered_png}" sig_zero.pdf
        parts="\$parts sig_zero.pdf"
        echo "[compile_master_report] Added zero-filtered signature plot"
    else
        echo "[compile_master_report] SKIP: zero_filtered signature PNG not available"
    fi
    if [ -s "${cosine_filtered_png}" ]; then
        Rscript /usr/local/bin/png_to_pdf.R "${cosine_filtered_png}" sig_cosine.pdf
        parts="\$parts sig_cosine.pdf"
        echo "[compile_master_report] Added cosine-filtered signature plot"
    else
        echo "[compile_master_report] SKIP: cosine_filtered signature PNG not available"
    fi

    # ── SECTION 4: Germline Prevalence Distributions ─────────────────────────
    add_section "Germline Prevalence Distributions" "sec04_prevalence.pdf"
    if [ -s "${prevalence_plot_png}" ]; then
        Rscript /usr/local/bin/png_to_pdf.R "${prevalence_plot_png}" prevalence.pdf
        parts="\$parts prevalence.pdf"
        echo "[compile_master_report] Added germline prevalence distributions plot"
    else
        echo "[compile_master_report] SKIP: germline prevalence plot not available"
    fi

    # ── SECTION 5: ADO Germline Summary ──────────────────────────────────────
    add_section "ADO Germline Summary" "sec05_ado_summary.pdf"
    if [ -s "${ado_summary_png}" ]; then
        Rscript /usr/local/bin/png_to_pdf.R "${ado_summary_png}" ado_summary.pdf
        parts="\$parts ado_summary.pdf"
        echo "[compile_master_report] Added ADO germline summary plot"
    else
        echo "[compile_master_report] SKIP: ADO summary plot not available"
    fi

    # ── SECTION 6: Phylogenetic Tree Topology Comparison ─────────────────────
    add_section "Phylogenetic Tree Topology Comparison" "sec06_tree_comparison.pdf"
    add_pdf "${tree_comparison_pdf}"

    # ── SECTION 7: Variant Filter Provenance Report ───────────────────────────
    add_section "Variant Filter Provenance Report" "sec07_filter_report.pdf"
    add_pdf "${combined_report_pdf}"

    # ── Assemble ──────────────────────────────────────────────────────────────
    echo "[compile_master_report] Parts: \$parts"
    pdfunite \$parts master_report_${group}.pdf
    echo "[compile_master_report] Done: master_report_${group}.pdf"
    """
}
