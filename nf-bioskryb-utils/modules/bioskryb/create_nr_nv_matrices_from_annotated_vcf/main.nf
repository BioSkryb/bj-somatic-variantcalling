nextflow.enable.dsl=2
params.timestamp    = ""
params.hq_min_pct   = 70   // % of samples that must have NR >= hq_min_depth for pileup_hq_depth scheme
params.hq_min_depth = 2    // minimum SMPL_PILEUP_NUM_FRAGMENTS_HQ_POSITION to count as "covered"

// Collect all per-sample annotated VCFs for a group and build cross-sample NR/NV matrices
// for 9 filtering schemes (unfiltered + 8 filtered). Each scheme outputs:
//   NR_annotated_vcf_${group}_<label>.tsv — rows=variants, cols=samples, values=SMPL_PILEUP_NUM_FRAGMENTS_HQ_POSITION
//   NV_annotated_vcf_${group}_<label>.tsv — rows=variants, cols=samples, values=HQ_MQ_BQ_F + HQ_MQ_BQ_R
// NR/NV values are extracted once (unfiltered); filtered schemes only differ in which variant rows are included.
// The union of variants passing the filter across all samples is used as the row index for each scheme.
process CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), val(sample_names), path(vcfs), path(tbis)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group),
          path("NR_annotated_vcf_${group}_*.tsv"),
          path("NV_annotated_vcf_${group}_*.tsv"),
          emit: nr_nv_matrices

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    first_vcf=\$(ls *_somatic_annotated.vcf.gz | sort | head -1)
    ls *_somatic_annotated.vcf.gz | sort | sed 's/_somatic_annotated\\.vcf\\.gz//' > sample_list.txt

    echo "Samples: \$(wc -l < sample_list.txt)"

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 1: Extract full per-sample NR and NV vectors (one value per variant)
    # ─────────────────────────────────────────────────────────────────────────
    bcftools query -f '%CHROM\\_%POS\\_%REF\\_%ALT\\n' "\${first_vcf}" > variant_ids_all.txt
    echo "Total variants: \$(wc -l < variant_ids_all.txt)"

    while IFS= read -r sample; do
        # NR: total HQ fragments at position (depth denominator)
        bcftools query \\
            -f '%INFO/SMPL_PILEUP_NUM_FRAGMENTS_HQ_POSITION\\n' \\
            "\${sample}_somatic_annotated.vcf.gz" \\
            | awk '{printf "%d\\n", (\$0+0)}' > "\${sample}_nr.txt"

        # NV: HQ ALT-supporting fragments (forward + reverse)
        bcftools query \\
            -f '%INFO/SMPL_PILEUP_NUM_FRAGMENTS_HQ_MQ_BQ_F\\t%INFO/SMPL_PILEUP_NUM_FRAGMENTS_HQ_MQ_BQ_R\\n' \\
            "\${sample}_somatic_annotated.vcf.gz" \\
            | awk -F'\\t' '{printf "%d\\n", (\$1+0)+(\$2+0)}' > "\${sample}_nv.txt"
    done < sample_list.txt

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
    }

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 2: Build unfiltered matrices (scheme 1 — baseline)
    # ─────────────────────────────────────────────────────────────────────────
    build_matrix variant_ids_all.txt unfiltered

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 3: Cohort-level filtered schemes (2–6)
    # All cohort-level INFO fields are identical across sample VCFs for a given
    # variant, so the first VCF is sufficient to derive the union.
    # ─────────────────────────────────────────────────────────────────────────

    # Scheme 2: Pileup_Verdict="Pass"
    bcftools query -i 'INFO/Pileup_Verdict="Pass"' \\
        -f '%CHROM\\_%POS\\_%REF\\_%ALT\\n' "\${first_vcf}" | sort -u > union_pileup.txt
    subset_matrix union_pileup.txt pileup

    # Scheme 3: Pileup_Verdict="Pass" AND ≥hq_min_pct% of samples have NR≥hq_min_depth
    # Thresholds come from params.hq_min_pct / params.hq_min_depth (default: 70 / 2).
    # Override in nextflow.config or with --hq_min_pct / --hq_min_depth on the command line.
    # Uses per-sample _nr.txt files already built in STEP 1.
    hq_min_pct=${params.hq_min_pct}
    hq_min_depth=${params.hq_min_depth}
    n_samples=\$(wc -l < sample_list.txt)

    # Paste all per-sample NR vectors alongside variant IDs into one wide table.
    paste variant_ids_all.txt \$(while IFS= read -r s; do echo "\${s}_nr.txt"; done < sample_list.txt) \\
        > all_nr_table.txt

    # Keep pileup variants where ≥hq_min_pct% of samples satisfy NR≥hq_min_depth.
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
    bcftools query -i 'INFO/Pileup_Verdict="Pass"' \\
        -f '%CHROM\\_%POS\\_%REF\\_%ALT\\t%INFO/SEQUOIA_SecondPass_Germline_qval_log10\\t%INFO/SEQUOIA_SecondPass_Rho\\n' \\
        "\${first_vcf}" \\
        | awk -F'\\t' '\$2 != "." && \$3 != "." && (\$2+0) < -5 && (\$3+0) > 0.2 {print \$1}' \\
        | sort -u > union_pileup_sequoia.txt
    subset_matrix union_pileup_sequoia.txt pileup_sequoia

    # Schemes 5–7: Pileup_Verdict_PassCount > N
    # Pileup_Verdict_PassCount is Type=String in the VCF header, so bcftools
    # arithmetic operators fail. Extract the field with bcftools query and
    # perform the numeric comparison in awk instead.
    for thresh in 1 2 3; do
        bcftools query -f '%CHROM\\_%POS\\_%REF\\_%ALT\\t%INFO/Pileup_Verdict_PassCount\\n' "\${first_vcf}" \\
            | awk -F'\\t' -v t="\${thresh}" '(\$2+0) > t {print \$1}' | sort -u > "union_passcount_gt\${thresh}.txt"
        subset_matrix "union_passcount_gt\${thresh}.txt" "pileup_passcount_gt\${thresh}"
    done

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 4: Per-sample filtered schemes (8–9)
    # SMPL_PILEUP_PropClipped_Filter varies per sample VCF — each sample is
    # queried independently and results are unioned.
    # ─────────────────────────────────────────────────────────────────────────

    # Scheme 8: SMPL_PILEUP_PropClipped_Filter="Pass" (union across samples)
    while IFS= read -r sample; do
        bcftools query -i 'INFO/SMPL_PILEUP_PropClipped_Filter="Pass"' \\
            -f '%CHROM\\_%POS\\_%REF\\_%ALT\\n' "\${sample}_somatic_annotated.vcf.gz"
    done < sample_list.txt | sort -u > union_smpl_propcl.txt
    subset_matrix union_smpl_propcl.txt smpl_propcl

    # Scheme 9: SMPL_PILEUP_PropClipped_Filter="Pass" AND SEQUOIA second-pass germline filters
    # SEQUOIA_SecondPass_* fields are Type=String — use awk for numeric comparisons.
    while IFS= read -r sample; do
        bcftools query -i 'INFO/SMPL_PILEUP_PropClipped_Filter="Pass"' \\
            -f '%CHROM\\_%POS\\_%REF\\_%ALT\\t%INFO/SEQUOIA_SecondPass_Germline_qval_log10\\t%INFO/SEQUOIA_SecondPass_Rho\\n' \\
            "\${sample}_somatic_annotated.vcf.gz"
    done < sample_list.txt \\
        | awk -F'\\t' '\$2 != "." && \$3 != "." && (\$2+0) < -5 && (\$3+0) > 0.2 {print \$1}' \\
        | sort -u > union_smpl_propcl_sequoia.txt
    subset_matrix union_smpl_propcl_sequoia.txt smpl_propcl_sequoia

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 5: Summary table — variant counts per scheme
    # ─────────────────────────────────────────────────────────────────────────
    printf "\\n%-40s %s\\n" "Matrix" "Variants"
    printf "%-40s %s\\n" "----------------------------------------" "--------"
    for label in unfiltered pileup pileup_hq_depth pileup_sequoia \\
                 pileup_passcount_gt1 pileup_passcount_gt2 pileup_passcount_gt3 \\
                 smpl_propcl smpl_propcl_sequoia; do
        nrows=\$(tail -n +2 NR_annotated_vcf_${group}_\${label}.tsv | wc -l)
        printf "%-40s %d\\n" "\${label}" "\${nrows}"
    done
    """
}
