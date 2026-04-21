nextflow.enable.dsl=2
params.timestamp = ""

// Extracts a long-format prevalence TSV from the IDENTIFY_GERMLINE_FROM_STATS
// annotated VCF for three filter categories:
//
//   GERMLINE_FROM_STATS   — GT-based non-REF prevalence for all Stage-1
//                           germline candidates (Binom q-value + Rho filter).
//                           Computed from VCF genotypes; covers the FULL
//                           distribution across the prevalence threshold.
//
//   VEP_AF_filter         — non-REF prevalence (stored INFO field) for all
//                           variants flagged by the VEP AF filter.
//
//   Bulk_Fail             — non-REF prevalence (stored INFO field) for all
//                           variants that failed the bulk filter.
//
// Output TSV columns: Variant  Filter  Prevalence_proportion
//
// The TSV is consumed by PLOT_GERMLINE_PREVALENCE_DISTRIBUTIONS.

process EXTRACT_GERMLINE_PREVALENCE_TABLE {
    tag "${group}"
    container 'quay.io/biocontainers/bcftools:1.14--h88f3f91_0'
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), path(vcf), path(tbi)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), path("germline_prevalence_long_${group}.tsv"), emit: prevalence_table

    script:
    """
    set -euo pipefail

    n_smp=\$(bcftools query -l ${vcf} | wc -l)
    echo "[EXTRACT_GERMLINE_PREVALENCE_TABLE] Samples: \${n_smp}"

    printf 'Variant\\tFilter\\tPrevalence_proportion\\n' \
        > germline_prevalence_long_${group}.tsv

    # ── 1. GERMLINE_FROM_STATS=Yes: GT-based prevalence (full distribution) ──
    echo "[EXTRACT_GERMLINE_PREVALENCE_TABLE] Querying GERMLINE_FROM_STATS variants..."
    bcftools query \
        -i 'INFO/GERMLINE_FROM_STATS="Yes"' \
        -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' \
        ${vcf} \
    | awk -v n_smp="\${n_smp}" -F'\\t' '{
        nonref = 0
        for (i = 5; i <= NF; i++) {
            gt = \$i; gsub(/[|]/, "/", gt)
            if (gt != "0/0" && gt != "./." && gt != ".") nonref++
        }
        printf "%s_%s_%s_%s\\tGERMLINE_FROM_STATS\\t%.6f\\n", \$1, \$2, \$3, \$4, nonref / n_smp
    }' >> germline_prevalence_long_${group}.tsv

    n_germ=\$(tail -n+2 germline_prevalence_long_${group}.tsv | grep -c GERMLINE || true)
    echo "[EXTRACT_GERMLINE_PREVALENCE_TABLE] GERMLINE_FROM_STATS rows: \${n_germ}"

    # ── 2. VEP_AF_filter: stored prevalence field ────────────────────────────
    echo "[EXTRACT_GERMLINE_PREVALENCE_TABLE] Querying VEP AF_filter variants..."
    bcftools query \
        -i 'INFO/VEP_AF_FILTER="AF_filter"' \
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/VEP_AF_FILTER_Prevalence_proportion\\n' \
        ${vcf} \
    | awk -F'\\t' '\$5 != "." { printf "%s_%s_%s_%s\\tVEP_AF_filter\\t%s\\n", \$1, \$2, \$3, \$4, \$5 }' \
    >> germline_prevalence_long_${group}.tsv

    n_vep=\$(tail -n+2 germline_prevalence_long_${group}.tsv | grep -c VEP || true)
    echo "[EXTRACT_GERMLINE_PREVALENCE_TABLE] VEP_AF_filter rows: \${n_vep}"

    # ── 3. Bulk_Fail: stored prevalence field ────────────────────────────────
    echo "[EXTRACT_GERMLINE_PREVALENCE_TABLE] Querying Bulk_Fail variants..."
    bcftools query \
        -i 'INFO/RemainingAfterBulk="Fail"' \
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/BULK_FAIL_Prevalence_proportion\\n' \
        ${vcf} \
    | awk -F'\\t' '\$5 != "." { printf "%s_%s_%s_%s\\tBulk_Fail\\t%s\\n", \$1, \$2, \$3, \$4, \$5 }' \
    >> germline_prevalence_long_${group}.tsv

    n_bulk=\$(tail -n+2 germline_prevalence_long_${group}.tsv | grep -c Bulk || true)
    echo "[EXTRACT_GERMLINE_PREVALENCE_TABLE] Bulk_Fail rows: \${n_bulk}"

    echo "[EXTRACT_GERMLINE_PREVALENCE_TABLE] Total rows (excl. header): \$(tail -n+2 germline_prevalence_long_${group}.tsv | wc -l)"
    """
}
