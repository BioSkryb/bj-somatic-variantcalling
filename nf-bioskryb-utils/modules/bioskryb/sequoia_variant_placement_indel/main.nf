nextflow.enable.dsl=2
params.timestamp = ""

// Place INDEL variants from an NR/NV matrix pair onto a pre-built phylogenetic
// tree using rscript_variant_placement.R (treemut::assign_to_tree).
//
// Default inputs (to be wired in the calling workflow):
//   tree_file  — pileup tree from SEQUOIA_PHYLOGENY_INDEL
//                e.g. output_indel_pileup/<group>_pileup_indel_for_MPBoot.fa.treefile
//   nr_matrix  — NR_annotated_vcf_<group>_unfiltered.tsv  (from CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF)
//   nv_matrix  — NV_annotated_vcf_<group>_unfiltered.tsv
//
// The R script handles the INDEL row subsetting internally (--variant_type indel).
// Gender is passed explicitly from params.gender (not inferred from coverage).
// The bash pre-checks mirror SEQUOIA_PHYLOGENY_INDEL: verify tree is non-empty,
// then require >=1 INDEL row and >=1 chrX/Y INDEL row before launching R.
process SEQUOIA_VARIANT_PLACEMENT_INDEL {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), path(tree_file), path(nr_matrix), path(nv_matrix)
    val(gender)
    val(vaf_absent)
    val(vaf_present)
    val(tree_mut_pval)
    val(keep_ancestral)
    val(create_multi_tree)
    val(genotype_conv_prob)
    val(min_pval_for_true_somatic)
    val(min_variant_reads_shared)
    val(min_vaf_shared)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), path("output_indel_placement_*"), emit: placement_outputs

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Sentinel so the output glob always resolves even when the run is skipped
    touch output_indel_placement_no_results

    # ── Guard: tree file must be non-empty ────────────────────────────────────
    tree_bytes=\$(wc -c < "${tree_file}" || echo 0)
    if [ "\${tree_bytes}" -eq 0 ]; then
        echo "[indel_placement] Tree file is empty (0 bytes) — skipping"
        exit 0
    fi

    # ── Count INDEL rows in the input NR matrix ───────────────────────────────
    n_indel=\$(awk 'NR>1 { n=split(\$1,a,"_"); ref=a[n-1]; alt=a[n];
                            if (length(ref)>1 || length(alt)>1) c++ }
                    END  { print c+0 }' "${nr_matrix}")
    echo "[indel_placement] INDEL variants in matrix: \${n_indel}"

    if [ "\${n_indel}" -eq 0 ]; then
        echo "[indel_placement] No INDELs — skipping"
        exit 0
    fi

    # ── Require at least one chrX/Y INDEL (needed for gender-aware VAF thresholds)
    n_sex=\$(awk 'NR>1 { n=split(\$1,a,"_"); ref=a[n-1]; alt=a[n];
                          if ((length(ref)>1 || length(alt)>1) &&
                              (\$1 ~ /^chrX_/ || \$1 ~ /^chrY_/)) c++ }
                  END  { print c+0 }' "${nr_matrix}")
    echo "[indel_placement] chrX/Y INDEL variants: \${n_sex}"

    if [ "\${n_sex}" -eq 0 ]; then
        echo "[indel_placement] No chrX/chrY INDELs — skipping"
        exit 0
    fi

    mkdir -p output_indel_placement_${group}

    Rscript /usr/local/bin/rscript_variant_placement.R \\
        --donor_id                  "${group}" \\
        --input_nr                  "${nr_matrix}" \\
        --input_nv                  "${nv_matrix}" \\
        --input_tree                "${tree_file}" \\
        --output_dir                "output_indel_placement_${group}/" \\
        --variant_type              "indel" \\
        --gender                    "${gender}" \\
        --vaf_absent                ${vaf_absent} \\
        --vaf_present               ${vaf_present} \\
        --tree_mut_pval             ${tree_mut_pval} \\
        --keep_ancestral            ${keep_ancestral} \\
        --create_multi_tree         ${create_multi_tree} \\
        --genotype_conv_prob        ${genotype_conv_prob} \\
        --min_pval_for_true_somatic ${min_pval_for_true_somatic} \\
        --min_variant_reads_shared  ${min_variant_reads_shared} \\
        --min_vaf_shared            ${min_vaf_shared}
    """
}
