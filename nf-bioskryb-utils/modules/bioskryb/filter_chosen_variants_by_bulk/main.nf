nextflow.enable.dsl=2
params.timestamp = ""

// Emit a single empty file when bulk VCF is not provided (so downstream filters/SEQUOIA get a valid channel).
process CREATE_EMPTY_BULK_VARIANTS {
    tag "empty_bulk"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: false

    input:
    val(publish_dir)
    val(enable_publish)

    output:
    path("empty_bulk_variants.txt"), emit: bulk_variants

    script:
    """
    touch empty_bulk_variants.txt
    """
}

// Remove from chosen_variants any variant that appears in the bulk variants list.
// bulk_variants_file may be empty (e.g. when bulk is not used); then all chosen variants pass through.
process FILTER_CHOSEN_VARIANTS_BY_BULK {
    tag "${group}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    tuple val(group), path(chosen_variants), path(bulk_variants)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(group), path("chosen_variants_filtered_bulk_${group}.txt"), emit: chosen_variants

    script:
    // When bulk file is empty, pass through all chosen variants. Otherwise awk: skip header in bulk; exclude chosen variants that appear in bulk.
    // (When first file is empty, awk's NR==FNR is true for every line of the second file, so all would be wrongly excluded.)
    """
    if [ ! -s "${bulk_variants}" ]; then
      cp ${chosen_variants} chosen_variants_filtered_bulk_${group}.txt
    else
      awk 'NR==FNR { if (FNR==1 && (\$0 ~ /^#/ || \$0 ~ /^CHROM/)) next; bulk[\$0]=1; next } !(\$0 in bulk)' ${bulk_variants} ${chosen_variants} > chosen_variants_filtered_bulk_${group}.txt
    fi
    """
}
