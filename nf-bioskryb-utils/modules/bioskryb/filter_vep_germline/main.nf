nextflow.enable.dsl=2
params.timestamp = ""

// Filter VEP-annotated VCF: optionally exclude variants only in dbSNP (Existing_variation);
// always apply AF/MAX_AF <= max_af. When filter_by_existing_variation is true: keep only if no IDs, or ≥1 non-dbSNP ID.
// Output: CHROM_POS_REF_ALT (one per line) to chosen_variants_postgermlinefilter_${group}.txt
process FILTER_VEP_GERMLINE {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), path(vep_vcf), path(vep_vcf_tbi)
    val(max_af)
    val(filter_by_existing_variation)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), path("chosen_variants_postgermlinefilter_${group}.txt"), emit: chosen_variants

    script:
    def vcf = vep_vcf[0]
    def filter_ev = filter_by_existing_variation ? 1 : 0
    // bcftools -i/-e does not see CSQ-extracted AF when using -f (format output), so do AF + optional Existing_variation in awk
    // -s worst: one row per variant (worst consequence). Cols: CHROM,POS,REF,ALT,Existing_variation,AF,MAX_AF
    """
    echo -e "Filtering VEP VCF: filter_by_existing_variation=${filter_by_existing_variation}; exclude AF, MAX_AF > ${max_af} ..."

    bcftools +split-vep ${vcf} -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%Existing_variation\\t%AF\\t%MAX_AF\\n' -d -s worst | \\
    awk -v max_af=${max_af} -v filter_ev=${filter_ev} -v OFS="_" 'BEGIN{FS="\\t"}
      NF>=7 {
        keep_ev=1;
        if (filter_ev) {
          ev=\$5; keep_ev=0;
          if (ev=="" || ev==".") keep_ev=1;
          else {
            gsub(/&/, ",", ev); n=split(ev,a,",");
            for(i=1;i<=n;i++) if (a[i]!="" && a[i]!~/^rs[0-9]+\$/) { keep_ev=1; break }
          }
        }
        if (keep_ev && (\$6=="" || \$6=="." || \$6+0<=max_af) && (\$7=="" || \$7=="." || \$7+0<=max_af)) {
          n=split(\$4,a,","); for(i=1;i<=n;i++) print \$1,\$2,\$3,a[i]
        }
      }' | sort -u > chosen_variants_postgermlinefilter_${group}.txt

    echo -e "Wrote \$(wc -l < chosen_variants_postgermlinefilter_${group}.txt) variants to chosen_variants_postgermlinefilter_${group}.txt"
    """
}
