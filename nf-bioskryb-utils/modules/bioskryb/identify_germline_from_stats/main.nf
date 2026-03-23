nextflow.enable.dsl=2
params.timestamp = ""

// Identify high-confidence germline variants from statistical evidence and annotate
// the merged group VCF with twelve new INFO fields in a single bcftools annotate pass:
//
//   HIGH_CONFIDENCE_GERMLINE_FROM_STATS
//       Stage 1 — annot_table: Binom_Germline_qval_log10 > -1 AND (Binom_Rho < 0.1 OR Binom_Rho = NA)
//       Stage 2 — group VCF genotypes: > 80 % of samples carry a het call (REF + any ALT, ploidy-agnostic)
//       Values: "Yes" | "No"
//
//   VEP_AF_FILTER
//       Variants in vep_filter_provenance with FILTER_REASON == "AF_filter" (exact match).
//       Values: "AF_filter" | "Not_filtered"
//
//   VEP_GERMLINE_HIGH_CONFIDENCE
//       Subset of VEP_AF_FILTER variants where >= 80 % of cells carry a het (0/1) call.
//       Values: "Yes" | "No"
//
//   RemainingAfterBulk
//       Pass/Fail status from bulk_filter_provenance (FILTER_CHOSEN_VARIANTS_BY_BULK output).
//       Variants absent from the provenance (not chosen) receive "Not_evaluated".
//       When no bulk VCF was used (empty provenance), all variants receive "Pass".
//
//   RemainingAfterBulk_HIGH_CONFIDENCE
//       Subset of RemainingAfterBulk=Fail variants where >= 80 % of cells carry a het call
//       (REF + any ALT, ploidy-agnostic).
//       Values: "Yes" | "No"
//
// Additional outputs:
//   germline_stats_prevalence_${group}.tsv — per-variant non-REF prevalence for HIGH_CONFIDENCE_GERMLINE_FROM_STATS=Yes variants
//   af_filter_prevalence_${group}.tsv      — per-variant non-REF prevalence for VEP AF_filter variants
//   gt_vep.txt                             — AF_filter variant IDs passing the 80 % het threshold
//   bulk_fail_prevalence_${group}.tsv      — per-variant non-REF prevalence for RemainingAfterBulk=Fail variants
//   gt_bulk.txt                            — Bulk Fail variant IDs passing the 80 % het threshold
process IDENTIFY_GERMLINE_FROM_STATS {
    tag "${group}"
    container 'quay.io/biocontainers/bcftools:1.14--h88f3f91_0'
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), path(group_vcf), path(annot_table), path(vep_filter_provenance), path(bulk_filter_provenance)
    val(germline_prev_pct)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group),
          path("${group}_hcgermline_annotated.vcf.gz"),
          path("${group}_hcgermline_annotated.vcf.gz.tbi"),
          emit: annotated_vcf

    script:
    def vcf = group_vcf instanceof List ? group_vcf[0] : group_vcf
    """
    set -euo pipefail

    # ── Stage 1: statistical filter on annot_table ────────────────────────────
    awk -F'\\t' '
    function qval_passes(v) {
        if (v == "NA" || v == "NaN") return 0
        if (v == "-Inf")             return 0
        if (v == "Inf")              return 1
        return v + 0 > -1
    }
    function rho_passes(v) {
        if (v == "NA" || v == "NaN") return 1
        if (v == "Inf" || v == "-Inf") return 0
        return v + 0 < 0.1
    }
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            if (\$i == "Binom_Germline_qval_log10") qcol = i
            if (\$i == "Binom_Rho")                 rcol = i
        }
        next
    }
    qcol && rcol && qval_passes(\$qcol) && rho_passes(\$rcol) { print \$1 }
    ' ${annot_table} | sort -u > candidate_ids.txt

    echo "[IDENTIFY_GERMLINE_FROM_STATS] Stage 1 candidates: \$(wc -l < candidate_ids.txt)"

    # ── Stage 2: genotype prevalence filter (> 80 % het) ─────────────────────
    n_samples=\$(bcftools query -l ${vcf} | wc -l)
    echo "[IDENTIFY_GERMLINE_FROM_STATS] Samples in VCF: \${n_samples}"

    touch germline_confirmed_ids_${group}.txt
    if [ -s "candidate_ids.txt" ]; then
        awk -F'_' '{
            n = NF; pos = \$(n-2)
            chrom = \$1; for (j = 2; j <= n-3; j++) chrom = chrom "_" \$j
            print chrom "\\t" (pos - 1) "\\t" pos
        }' candidate_ids.txt | sort -k1,1 -k2,2n | uniq > candidate_regions.bed

        bcftools query \
            -R candidate_regions.bed \
            -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' \
            ${vcf} > candidate_gts.tsv

        awk -v n_smp="\${n_samples}" -v pct=${germline_prev_pct} -F'\\t' '
        NR == FNR {
            vid = \$0; n = split(vid, a, "_")
            alt = a[n]; ref = a[n-1]; pos = a[n-2]
            chrom = a[1]; for (j = 2; j <= n-3; j++) chrom = chrom "_" a[j]
            cands[chrom "_" pos "_" ref "_" alt] = vid
            next
        }
        {
            key = \$1 "_" \$2 "_" \$3 "_" \$4
            if (!(key in cands)) next
            het = 0
            for (i = 5; i <= NF; i++) {
                gt = \$i; gsub(/[|]/, "/", gt)
                n_al = split(gt, al, "/")
                has_ref = 0; has_alt = 0
                for (a = 1; a <= n_al; a++) {
                    if (al[a] == "0") has_ref = 1
                    else if (al[a] != ".") has_alt = 1
                }
                if (has_ref && has_alt) het++
            }
            if (het / n_smp * 100 > pct) print cands[key]
        }
        ' candidate_ids.txt candidate_gts.tsv | sort -u > germline_confirmed_ids_${group}.txt
    fi

    echo "[IDENTIFY_GERMLINE_FROM_STATS] HIGH_CONFIDENCE_GERMLINE_FROM_STATS=Yes: \$(wc -l < germline_confirmed_ids_${group}.txt)"

    # ── HIGH_CONFIDENCE_GERMLINE_FROM_STATS: prevalence table ────────────────
    # Prevalence = cells where GT is not 0/0 (any non-REF call, including het and hom-alt).
    # candidate_gts.tsv (built in Stage 2) covers all Stage 1 candidates — a superset
    # of confirmed IDs — so no extra bcftools query is needed here.
    printf 'Variant\\tPrevalence_count\\tPrevalence_proportion\\n' \
        > germline_stats_prevalence_${group}.tsv

    if [ -s "germline_confirmed_ids_${group}.txt" ]; then
        awk -v n_smp="\${n_samples}" -F'\\t' '
        NR == FNR { conf[\$0] = 1; next }
        {
            key = \$1 "_" \$2 "_" \$3 "_" \$4
            if (!(key in conf)) next
            nonref = 0
            for (i = 5; i <= NF; i++) {
                gt = \$i; gsub(/[|]/, "/", gt)
                if (gt != "0/0" && gt != "./." && gt != ".") nonref++
            }
            printf "%s\\t%d\\t%.4f\\n", key, nonref, nonref / n_smp
        }
        ' germline_confirmed_ids_${group}.txt candidate_gts.tsv \
            >> germline_stats_prevalence_${group}.tsv
    fi

    echo "[IDENTIFY_GERMLINE_FROM_STATS] germline_stats_prevalence rows: \$(tail -n+2 germline_stats_prevalence_${group}.tsv | wc -l)"

    # ── Parse VEP AF_filter variants ──────────────────────────────────────────
    if [ -s "${vep_filter_provenance}" ]; then
        awk -F'\\t' '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if (\$i == "CHROM")         cc = i
                if (\$i == "POS")           pc = i
                if (\$i == "REF")           rc = i
                if (\$i == "ALT")           ac = i
                if (\$i == "FILTER_REASON") fc = i
            }
            next
        }
        cc && fc && \$fc == "AF_filter" { print \$cc "_" \$pc "_" \$rc "_" \$ac }
        ' ${vep_filter_provenance} | sort -u > af_filter_ids.txt
    else
        touch af_filter_ids.txt
    fi
    echo "[IDENTIFY_GERMLINE_FROM_STATS] AF_filter variants from VEP: \$(wc -l < af_filter_ids.txt)"

    # ── AF_filter: prevalence table + VEP_GERMLINE_HIGH_CONFIDENCE (>= 80 % het) ──
    # Prevalence = cells where GT is not REF (0/0) and not missing (./.)
    # VEP_GERMLINE_HIGH_CONFIDENCE = Yes when >= 80 % of cells have 0/1 het call.
    touch gt_vep.txt
    printf 'Variant\\tPrevalence_count\\tPrevalence_proportion\\n' \
        > af_filter_prevalence_${group}.tsv

    if [ -s "af_filter_ids.txt" ]; then
        awk -F'_' '{
            n = NF; pos = \$(n-2)
            chrom = \$1; for (j = 2; j <= n-3; j++) chrom = chrom "_" \$j
            print chrom "\\t" (pos - 1) "\\t" pos
        }' af_filter_ids.txt | sort -k1,1 -k2,2n | uniq > af_filter_regions.bed

        bcftools query \
            -R af_filter_regions.bed \
            -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' \
            ${vcf} > af_filter_gts.tsv

        awk -v n_smp="\${n_samples}" -v pct=${germline_prev_pct} -F'\\t' '
        NR == FNR { af_ids[\$0] = 1; next }
        {
            key = \$1 "_" \$2 "_" \$3 "_" \$4
            if (!(key in af_ids)) next
            nonref = 0; het = 0
            for (i = 5; i <= NF; i++) {
                gt = \$i; gsub(/[|]/, "/", gt)
                if (gt != "0/0" && gt != "./." && gt != ".") nonref++
                n_al = split(gt, al, "/")
                has_ref = 0; has_alt = 0
                for (a = 1; a <= n_al; a++) {
                    if (al[a] == "0") has_ref = 1
                    else if (al[a] != ".") has_alt = 1
                }
                if (has_ref && has_alt) het++
            }
            printf "%s\\t%d\\t%.4f\\n", key, nonref, nonref / n_smp
            if (het / n_smp * 100 >= pct) print key > "gt_vep.txt"
        }
        ' af_filter_ids.txt af_filter_gts.tsv >> af_filter_prevalence_${group}.tsv
    fi

    echo "[IDENTIFY_GERMLINE_FROM_STATS] VEP_GERMLINE_HIGH_CONFIDENCE=Yes: \$(wc -l < gt_vep.txt)"

    # ── Parse bulk filter provenance ──────────────────────────────────────────
    # Format: Variant<TAB>RemainingAfterBulk  (header on line 1)
    # When no bulk VCF was used the file is empty → all variants receive "Pass".
    # Variants present in the VCF but absent from the provenance were never chosen
    # variants and receive "Not_evaluated".
    if [ -s "${bulk_filter_provenance}" ]; then
        awk -F'\\t' 'NR > 1 && NF == 2 { bulk_status[\$1] = \$2 }
        END { for (v in bulk_status) print v "\\t" bulk_status[v] }' \
            ${bulk_filter_provenance} > bulk_lookup.tsv
        bulk_mode="provenance"
    else
        touch bulk_lookup.tsv
        bulk_mode="no_bulk"
    fi
    echo "[IDENTIFY_GERMLINE_FROM_STATS] Bulk mode: \${bulk_mode} (\$(wc -l < bulk_lookup.tsv) entries)"

    # ── Bulk Fail: prevalence table + RemainingAfterBulk_HIGH_CONFIDENCE (>= 80 % het) ──
    # Prevalence = cells where GT is not REF (0/0) and not missing (./.)
    # RemainingAfterBulk_HIGH_CONFIDENCE = Yes when >= 80 % of cells have a het call.
    touch gt_bulk.txt
    printf 'Variant\\tPrevalence_count\\tPrevalence_proportion\\n' \
        > bulk_fail_prevalence_${group}.tsv

    if [ -s "bulk_lookup.tsv" ]; then
        awk -F'\\t' '\$2 == "Fail" { print \$1 }' bulk_lookup.tsv | sort -u > bulk_fail_ids.txt
    else
        touch bulk_fail_ids.txt
    fi
    echo "[IDENTIFY_GERMLINE_FROM_STATS] RemainingAfterBulk=Fail variants: \$(wc -l < bulk_fail_ids.txt)"

    if [ -s "bulk_fail_ids.txt" ]; then
        awk -F'_' '{
            n = NF; pos = \$(n-2)
            chrom = \$1; for (j = 2; j <= n-3; j++) chrom = chrom "_" \$j
            print chrom "\\t" (pos - 1) "\\t" pos
        }' bulk_fail_ids.txt | sort -k1,1 -k2,2n | uniq > bulk_fail_regions.bed

        bcftools query \
            -R bulk_fail_regions.bed \
            -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' \
            ${vcf} > bulk_fail_gts.tsv

        awk -v n_smp="\${n_samples}" -v pct=${germline_prev_pct} -F'\\t' '
        NR == FNR { bulk_ids[\$0] = 1; next }
        {
            key = \$1 "_" \$2 "_" \$3 "_" \$4
            if (!(key in bulk_ids)) next
            nonref = 0; het = 0
            for (i = 5; i <= NF; i++) {
                gt = \$i; gsub(/[|]/, "/", gt)
                if (gt != "0/0" && gt != "./." && gt != ".") nonref++
                n_al = split(gt, al, "/")
                has_ref = 0; has_alt = 0
                for (a = 1; a <= n_al; a++) {
                    if (al[a] == "0") has_ref = 1
                    else if (al[a] != ".") has_alt = 1
                }
                if (has_ref && has_alt) het++
            }
            printf "%s\\t%d\\t%.4f\\n", key, nonref, nonref / n_smp
            if (het / n_smp * 100 >= pct) print key > "gt_bulk.txt"
        }
        ' bulk_fail_ids.txt bulk_fail_gts.tsv >> bulk_fail_prevalence_${group}.tsv
    fi

    echo "[IDENTIFY_GERMLINE_FROM_STATS] RemainingAfterBulk_HIGH_CONFIDENCE=Yes: \$(wc -l < gt_bulk.txt)"

    # ── Build combined annotation TSV (all VCF variants, five fields) ────────
    # Columns: CHROM POS REF ALT HIGH_CONFIDENCE_GERMLINE_FROM_STATS VEP_AF_FILTER VEP_GERMLINE_HIGH_CONFIDENCE RemainingAfterBulk RemainingAfterBulk_HIGH_CONFIDENCE
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${vcf} | \
    awk -F'\\t' -v bulk_mode="\${bulk_mode}" '
    BEGIN {
        while ((getline line < "germline_confirmed_ids_${group}.txt") > 0) {
            vid = line; n = split(vid, a, "_")
            alt = a[n]; ref = a[n-1]; pos = a[n-2]
            chrom = a[1]; for (j = 2; j <= n-3; j++) chrom = chrom "_" a[j]
            confirmed[chrom "_" pos "_" ref "_" alt] = 1
        }
        while ((getline line < "af_filter_ids.txt") > 0)
            af_ids[line] = 1
        while ((getline line < "gt_vep.txt") > 0)
            gt_vep_ids[line] = 1
        while ((getline line < "gt_bulk.txt") > 0)
            gt_bulk_ids[line] = 1
        while ((getline line < "bulk_lookup.tsv") > 0) {
            n = split(line, a, "\\t")
            if (n == 2) bulk[a[1]] = a[2]
        }
        while ((getline line < "candidate_ids.txt") > 0) {
            vid = line; n = split(vid, a, "_")
            alt = a[n]; ref = a[n-1]; pos = a[n-2]
            chrom = a[1]; for (j = 2; j <= n-3; j++) chrom = chrom "_" a[j]
            stage1_ids[chrom "_" pos "_" ref "_" alt] = 1
        }
        while ((getline line < "germline_stats_prevalence_${group}.tsv") > 0) {
            n = split(line, a, "\\t")
            if (n == 3 && a[1] != "Variant") { gstat_c[a[1]] = a[2]; gstat_p[a[1]] = a[3] }
        }
        while ((getline line < "af_filter_prevalence_${group}.tsv") > 0) {
            n = split(line, a, "\\t")
            if (n == 3 && a[1] != "Variant") { vep_c[a[1]] = a[2]; vep_p[a[1]] = a[3] }
        }
        while ((getline line < "bulk_fail_prevalence_${group}.tsv") > 0) {
            n = split(line, a, "\\t")
            if (n == 3 && a[1] != "Variant") { bkf_c[a[1]] = a[2]; bkf_p[a[1]] = a[3] }
        }
    }
    {
        key = \$1 "_" \$2 "_" \$3 "_" \$4
        hcg     = (key in confirmed)   ? "Yes"       : "No"
        af      = (key in af_ids)      ? "AF_filter"  : "Not_filtered"
        vep_hc  = (key in gt_vep_ids)  ? "Yes"        : "No"
        if (bulk_mode == "no_bulk")      rab = "Pass"
        else if (key in bulk)            rab = bulk[key]
        else                             rab = "Not_evaluated"
        rab_hc  = (key in gt_bulk_ids) ? "Yes"        : "No"
        gfs     = (key in stage1_ids)  ? "Yes"        : "No"
        gs_cnt  = (key in gstat_c)     ? gstat_c[key] : "."
        gs_pct  = (key in gstat_c)     ? gstat_p[key] : "."
        vaf_cnt = (key in vep_c)       ? vep_c[key]   : "."
        vaf_pct = (key in vep_c)       ? vep_p[key]   : "."
        bk_cnt  = (key in bkf_c)       ? bkf_c[key]   : "."
        bk_pct  = (key in bkf_c)       ? bkf_p[key]   : "."
        print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" hcg "\\t" af "\\t" vep_hc "\\t" rab "\\t" rab_hc "\\t" gfs "\\t" gs_cnt "\\t" gs_pct "\\t" vaf_cnt "\\t" vaf_pct "\\t" bk_cnt "\\t" bk_pct
    }' | sort -k1,1 -k2,2n | bgzip -c > annot_combined.tsv.gz
    tabix -s1 -b2 -e2 annot_combined.tsv.gz

    # ── Write INFO header lines ───────────────────────────────────────────────
    printf '##INFO=<ID=HIGH_CONFIDENCE_GERMLINE_FROM_STATS,Number=1,Type=String,Description="Germline: Binom_Germline_qval_log10>-1 AND Binom_Rho<0.1 AND >%d%% samples het (0/1)">\n' \
        "${germline_prev_pct}" > combined_hdr.txt
    printf '##INFO=<ID=VEP_AF_FILTER,Number=1,Type=String,Description="AF_filter if variant removed by VEP germline filter (FILTER_REASON==AF_filter exact); Not_filtered otherwise">\n' \
        >> combined_hdr.txt
    printf '##INFO=<ID=VEP_GERMLINE_HIGH_CONFIDENCE,Number=1,Type=String,Description="Yes if VEP AF_filter variant with >=%d%% of cells carrying het (0/1) genotype">\n' \
        "${germline_prev_pct}" >> combined_hdr.txt
    printf '##INFO=<ID=RemainingAfterBulk,Number=1,Type=String,Description="Bulk filter status: Pass/Fail from bulk_filter_provenance; Not_evaluated if variant was not a chosen variant; Pass when no bulk VCF was used">\n' \
        >> combined_hdr.txt
    printf '##INFO=<ID=RemainingAfterBulk_HIGH_CONFIDENCE,Number=1,Type=String,Description="Yes if RemainingAfterBulk=Fail variant with >=%d%% of cells carrying a het call (REF + any ALT, ploidy-agnostic)">\n' \
        "${germline_prev_pct}" >> combined_hdr.txt
    printf '##INFO=<ID=GERMLINE_FROM_STATS,Number=1,Type=String,Description="Stage 1 germline: Binom_Germline_qval_log10>-1 AND (Binom_Rho<0.1 OR Binom_Rho=NA); no genotype prevalence filter">\n' \
        >> combined_hdr.txt
    printf '##INFO=<ID=GERMLINE_STATS_Prevalence_count,Number=1,Type=Integer,Description="Cells with non-REF genotype for HIGH_CONFIDENCE_GERMLINE_FROM_STATS=Yes variants; missing for all others">\n' \
        >> combined_hdr.txt
    printf '##INFO=<ID=GERMLINE_STATS_Prevalence_proportion,Number=1,Type=Float,Description="Proportion of cells with non-REF genotype for HIGH_CONFIDENCE_GERMLINE_FROM_STATS=Yes variants; missing for all others">\n' \
        >> combined_hdr.txt
    printf '##INFO=<ID=VEP_AF_FILTER_Prevalence_count,Number=1,Type=Integer,Description="Cells with non-REF genotype for VEP_AF_FILTER=AF_filter variants; missing for all others">\n' \
        >> combined_hdr.txt
    printf '##INFO=<ID=VEP_AF_FILTER_Prevalence_proportion,Number=1,Type=Float,Description="Proportion of cells with non-REF genotype for VEP_AF_FILTER=AF_filter variants; missing for all others">\n' \
        >> combined_hdr.txt
    printf '##INFO=<ID=BULK_FAIL_Prevalence_count,Number=1,Type=Integer,Description="Cells with non-REF genotype for RemainingAfterBulk=Fail variants; missing for all others">\n' \
        >> combined_hdr.txt
    printf '##INFO=<ID=BULK_FAIL_Prevalence_proportion,Number=1,Type=Float,Description="Proportion of cells with non-REF genotype for RemainingAfterBulk=Fail variants; missing for all others">\n' \
        >> combined_hdr.txt

    # ── Single bcftools annotate pass (twelve INFO fields) ───────────────────
    bcftools annotate \
        -a annot_combined.tsv.gz \
        -h combined_hdr.txt \
        -c CHROM,POS,REF,ALT,INFO/HIGH_CONFIDENCE_GERMLINE_FROM_STATS,INFO/VEP_AF_FILTER,INFO/VEP_GERMLINE_HIGH_CONFIDENCE,INFO/RemainingAfterBulk,INFO/RemainingAfterBulk_HIGH_CONFIDENCE,INFO/GERMLINE_FROM_STATS,INFO/GERMLINE_STATS_Prevalence_count,INFO/GERMLINE_STATS_Prevalence_proportion,INFO/VEP_AF_FILTER_Prevalence_count,INFO/VEP_AF_FILTER_Prevalence_proportion,INFO/BULK_FAIL_Prevalence_count,INFO/BULK_FAIL_Prevalence_proportion \
        -Oz -o ${group}_hcgermline_annotated.vcf.gz \
        ${vcf}
    bcftools index -t ${group}_hcgermline_annotated.vcf.gz

    echo "[IDENTIFY_GERMLINE_FROM_STATS] Done — final variant counts:"
    bcftools view -H ${group}_hcgermline_annotated.vcf.gz \
        | awk '
          \$8 ~ /HIGH_CONFIDENCE_GERMLINE_FROM_STATS=Yes/ { hcg++ }
          \$8 ~ /GERMLINE_FROM_STATS=Yes/                 { gfs++ }
          \$8 ~ /VEP_AF_FILTER=AF_filter/                 { af++ }
          \$8 ~ /VEP_GERMLINE_HIGH_CONFIDENCE=Yes/        { vhc++ }
          \$8 ~ /RemainingAfterBulk=Pass/                 { rbp++ }
          \$8 ~ /RemainingAfterBulk=Fail/                 { rbf++ }
          \$8 ~ /RemainingAfterBulk=Not_evaluated/        { rbn++ }
          \$8 ~ /RemainingAfterBulk_HIGH_CONFIDENCE=Yes/  { rab_hc++ }
          END {
            print "  HIGH_CONFIDENCE_GERMLINE_FROM_STATS=Yes    : " hcg+0
            print "  GERMLINE_FROM_STATS=Yes                    : " gfs+0
            print "  VEP_AF_FILTER=AF_filter                    : " af+0
            print "  VEP_GERMLINE_HIGH_CONFIDENCE=Yes           : " vhc+0
            print "  RemainingAfterBulk=Pass                    : " rbp+0
            print "  RemainingAfterBulk=Fail                    : " rbf+0
            print "  RemainingAfterBulk=Not_evaluated           : " rbn+0
            print "  RemainingAfterBulk_HIGH_CONFIDENCE=Yes     : " rab_hc+0
          }'
    """
}
