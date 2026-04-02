nextflow.enable.dsl=2
params.timestamp               = ""
params.hq_min_pct              = 70   // % of samples that must have NR >= hq_min_depth for pileup_hq_depth scheme
params.hq_min_depth            = 2    // minimum SMPL_PILEUP_NUM_FRAGMENTS_HQ_POSITION to count as "covered"
params.second_run_hq_binom     = -5   // SEQUOIA_SecondPass_Germline_qval_log10 threshold for pileup_sequoia
params.second_run_hq_betabinom = 0.2  // SEQUOIA_SecondPass_Rho threshold for pileup_sequoia
// Only four schemes are built: unfiltered, pileup, pileup_hq_depth, pileup_sequoia.
// Optional scheme gating (params.nr_nv_matrix_schemes + scheme_requested) is disabled; see commented blocks in script.
params.nr_nv_matrix_schemes = ""

// Assemble cross-sample NR/NV/GT matrices from pre-extracted per-sample vectors.
// Input: grouped per-sample NR/NV/GT txt files produced by EXTRACT_NR_NV_GT_FROM_ANNOTATED_VCF.
// Each scheme outputs:
//   NR_annotated_vcf_${group}_<label>.tsv — rows=variants, cols=samples, values=SMPL_PILEUP_NUM_FRAGMENTS_HQ_POSITION
//   NV_annotated_vcf_${group}_<label>.tsv — rows=variants, cols=samples, values=HQ_MQ_BQ_F + HQ_MQ_BQ_R
// NR/NV values are extracted once (unfiltered); filtered schemes only differ in which variant rows are included.
// The union of variants passing the filter across all samples is used as the row index for each scheme.
// Pileup_Verdict and SEQUOIA_SecondPass_* fields are read from the first annotated VCF (group-level INFO fields
// are identical across samples for a given variant).
process CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), val(sample_names),
          path(nr_files),
          path(nv_files),
          path(gt_files),
          path(variant_ids_files),
          path(vcfs),
          path(tbis)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group),
          path("NR_annotated_vcf_${group}_*.tsv"),
          path("NV_annotated_vcf_${group}_*.tsv"),
          emit: nr_nv_matrices
    tuple val(group),
          path("matrix_scheme_summary_${group}.tsv"),
          emit: matrix_scheme_summary
    tuple val(group),
          path("matrix_per_sample_summary_${group}.tsv"),
          emit: matrix_per_sample_summary

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 1: Build sample list and collect pre-extracted per-sample vectors
    # ─────────────────────────────────────────────────────────────────────────
    # Sample list — derived from the nr file names (<sample>_nr.txt)
    ls *_nr.txt | sort | sed 's/_nr\\.txt//' > sample_list.txt
    echo "Samples: \$(wc -l < sample_list.txt)"

    # Variant IDs — all *_variant_ids.txt files are identical; use any one
    cp \$(ls *_variant_ids.txt | head -1) variant_ids_all.txt
    echo "Total variants: \$(wc -l < variant_ids_all.txt)"

    # First annotated VCF — used for cohort-level INFO field queries (Pileup_Verdict, SEQUOIA_SecondPass_*)
    first_vcf=\$(ls *_somatic_annotated*.vcf.gz | sort | head -1)

    echo "NR/NV matrix schemes: unfiltered, pileup, pileup_hq_depth, pileup_sequoia"

    # ─────────────────────────────────────────────────────────────────────────
    # Helper: assemble NR and NV matrices from a variant ID list
    # Usage: build_matrix <variant_ids_file> <label>
    # ─────────────────────────────────────────────────────────────────────────
    build_matrix() {
        local ids=\$1 label=\$2
        {
            printf ""
            while IFS= read -r s; do printf "\\t%s" "\$s"; done < sample_list.txt
            printf "\\n"
            paste "\$ids" \$(while IFS= read -r s; do echo "\${s}_nr.txt"; done < sample_list.txt)
        } > NR_annotated_vcf_${group}_\${label}.tsv

        {
            printf ""
            while IFS= read -r s; do printf "\\t%s" "\$s"; done < sample_list.txt
            printf "\\n"
            paste "\$ids" \$(while IFS= read -r s; do echo "\${s}_nv.txt"; done < sample_list.txt)
        } > NV_annotated_vcf_${group}_\${label}.tsv
    }

    # Helper: build GT matrix from a variant ID list
    # Usage: build_gt_matrix <variant_ids_file> <label>
    build_gt_matrix() {
        local ids=\$1 label=\$2
        {
            printf ""
            while IFS= read -r s; do printf "\\t%s" "\$s"; done < sample_list.txt
            printf "\\n"
            paste "\$ids" \$(while IFS= read -r s; do echo "\${s}_gt.txt"; done < sample_list.txt)
        } > GT_annotated_vcf_${group}_\${label}.tsv
    }

    # Helper: subset unfiltered matrices to variant rows present in a union file
    # Usage: subset_matrix <union_file> <label>
    subset_matrix() {
        local union=\$1 label=\$2
        awk 'NR==FNR {ids[\$1]=1; next} FNR==1 || (\$1 in ids)' \\
            "\$union" NR_annotated_vcf_${group}_unfiltered.tsv \\
            > NR_annotated_vcf_${group}_\${label}.tsv
        awk 'NR==FNR {ids[\$1]=1; next} FNR==1 || (\$1 in ids)' \\
            "\$union" NV_annotated_vcf_${group}_unfiltered.tsv \\
            > NV_annotated_vcf_${group}_\${label}.tsv
        awk 'NR==FNR {ids[\$1]=1; next} FNR==1 || (\$1 in ids)' \\
            "\$union" GT_annotated_vcf_${group}_unfiltered.tsv \\
            > GT_annotated_vcf_${group}_\${label}.tsv
    }

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 2: Build unfiltered matrices (always — required by downstream processes)
    # ─────────────────────────────────────────────────────────────────────────
    build_matrix variant_ids_all.txt unfiltered
    build_gt_matrix variant_ids_all.txt unfiltered

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 3: Cohort-level filtered schemes (2–4)
    # All cohort-level INFO fields are identical across sample VCFs for a given
    # variant, so the first VCF is sufficient to derive the union.
    # union_pileup.txt and pileup matrices are always built (required downstream).
    # ─────────────────────────────────────────────────────────────────────────

    # Scheme 2: Pileup_Verdict="Pass" — always built
    bcftools query -i 'INFO/Pileup_Verdict="Pass"' \\
        -f '%CHROM\\_%POS\\_%REF\\_%ALT\\n' "\${first_vcf}" | sort -u > union_pileup.txt
    subset_matrix union_pileup.txt pileup

    # Scheme 3: Pileup_Verdict="Pass" AND >=hq_min_pct% of samples have NR>=hq_min_depth
    # Thresholds come from params.hq_min_pct / params.hq_min_depth (default: 70 / 2).
    hq_min_pct=${params.hq_min_pct}
    hq_min_depth=${params.hq_min_depth}
    n_samples=\$(wc -l < sample_list.txt)
    paste variant_ids_all.txt \$(while IFS= read -r s; do echo "\${s}_nr.txt"; done < sample_list.txt) \\
        > all_nr_table.txt
    awk -F'\\t' -v pct="\${hq_min_pct}" -v depth="\${hq_min_depth}" -v n="\${n_samples}" '
        NR==FNR { pileup_ids[\$1]=1; next }
        !(\$1 in pileup_ids) { next }
        {
            pass_count = 0
            for (i=2; i<=NF; i++) if (\$i+0 >= depth) pass_count++
            if (pass_count / n * 100 >= pct) print \$1
        }
    ' union_pileup.txt all_nr_table.txt | sort -u > union_pileup_hq_depth.txt
    subset_matrix union_pileup_hq_depth.txt pileup_hq_depth

    # Scheme 4: Pileup_Verdict="Pass" AND SEQUOIA second-pass germline filters
    # SEQUOIA_SecondPass_* fields are Type=String — use awk for numeric comparisons.
    # Thresholds: params.second_run_hq_binom (default -5) and params.second_run_hq_betabinom (default 0.2)
    bcftools query -i 'INFO/Pileup_Verdict="Pass"' \\
        -f '%CHROM\\_%POS\\_%REF\\_%ALT\\t%INFO/SEQUOIA_SecondPass_Germline_qval_log10\\t%INFO/SEQUOIA_SecondPass_Rho\\n' \\
        "\${first_vcf}" \\
        | awk -F'\\t' '\$2 != "." && \$3 != "." && (\$2+0) < ${params.second_run_hq_binom} && (\$3+0) > ${params.second_run_hq_betabinom} {print \$1}' \\
        | sort -u > union_pileup_sequoia.txt
    subset_matrix union_pileup_sequoia.txt pileup_sequoia

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 5: Summary table — variant counts per scheme (only built schemes)
    # ─────────────────────────────────────────────────────────────────────────
    printf "\\n%-40s %s\\n" "Matrix" "Variants"
    printf "%-40s %s\\n" "----------------------------------------" "--------"
    for nr_file in NR_annotated_vcf_${group}_*.tsv; do
        label=\$(basename "\${nr_file}" | sed 's/NR_annotated_vcf_${group}_//' | sed 's/\\.tsv\$//')
        nrows=\$(tail -n +2 "\${nr_file}" | wc -l)
        printf "%-40s %d\\n" "\${label}" "\${nrows}"
    done

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 6: SNV vs INDEL counts per scheme
    # ─────────────────────────────────────────────────────────────────────────
    count_snv_indel() {
        local ids_file=\$1
        if [ ! -s "\${ids_file}" ]; then
            printf "0\\t0\\n"; return
        fi
        awk '{
            n = split(\$0, a, "_")
            ref = a[n-1]; alt = a[n]
            if (length(ref) == 1 && length(alt) == 1) snv++
            else indel++
        } END { printf "%d\\t%d\\n", snv+0, indel+0 }' "\${ids_file}"
    }

    printf "scheme\\tNumberOfSNVs\\tNumberOfIndels\\n" > matrix_scheme_summary_${group}.tsv

    for label in unfiltered pileup pileup_hq_depth pileup_sequoia; do
        case "\${label}" in
            unfiltered)      ids=variant_ids_all.txt ;;
            pileup)          ids=union_pileup.txt ;;
            pileup_hq_depth) ids=union_pileup_hq_depth.txt ;;
            pileup_sequoia)  ids=union_pileup_sequoia.txt ;;
        esac
        counts=\$(count_snv_indel "\${ids}")
        printf "%s\\t%s\\n" "\${label}" "\${counts}" >> matrix_scheme_summary_${group}.tsv
    done

    echo "[matrix_scheme_summary] Done:"
    cat matrix_scheme_summary_${group}.tsv

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 7: Per-sample SNV / INDEL counts per scheme using GT matrices
    # ─────────────────────────────────────────────────────────────────────────
    printf "scheme\\tsample\\tNumberOfSNVs\\tNumberOfIndels\\n" > matrix_per_sample_summary_${group}.tsv

    for gt_file in GT_annotated_vcf_${group}_*.tsv; do
        label=\$(basename "\${gt_file}" | sed 's/GT_annotated_vcf_${group}_//' | sed 's/\\.tsv\$//')
        awk -F'\\t' -v scheme="\${label}" '
        NR==1 {
            for (i=2; i<=NF; i++) samples[i]=\$i
            ncols=NF; next
        }
        {
            n = split(\$1, a, "_")
            ref = a[n-1]; alt = a[n]
            is_snv = (length(ref)==1 && length(alt)==1)
            for (i=2; i<=ncols; i++) {
                gt = \$i
                if (gt != "0/0" && gt != "0|0" && gt != "./." && gt != ".|." && gt != ".") {
                    if (is_snv) snv[i]++
                    else        indel[i]++
                }
            }
        }
        END {
            for (i=2; i<=ncols; i++)
                printf "%s\\t%s\\t%d\\t%d\\n", scheme, samples[i], snv[i]+0, indel[i]+0
        }' "\${gt_file}" >> matrix_per_sample_summary_${group}.tsv
    done

    echo "[matrix_per_sample_summary] Done: \$(tail -n +2 matrix_per_sample_summary_${group}.tsv | wc -l) rows"
    """
}
