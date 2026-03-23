nextflow.enable.dsl=2
params.timestamp = ""

// Variant of CONCAT_SUMMARY_ADO_INTERVALS that accepts a provenance label
// ("stats", "vep", or "bulk") and names all outputs with that label, so
// per-provenance summaries are kept separate in the publish directory.
//
// Inputs are the res_ADO_* TSV files collected from SUMMARIZE_ADO_INTERVALS
// after they have been split by provenance (see pipeline wiring).
//
// Output files:
//   merged_ADO_${label}.tsv        — combined allele-frequency distribution
//   ADO_plot_summary_${label}.png  — summary plot
//   merged_ADO_summary_${label}.tsv — per-sample ADO percentage

process CONCAT_SUMMARY_ADO_INTERVALS_LABELED {
    tag "${label}"
    container '597246834581.dkr.ecr.us-east-1.amazonaws.com/miscellaneous:r_ado_1.2'
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled: "$enable_publish"

    input:
    path(df_files)
    val(label)
    val(publish_dir)
    val(enable_publish)

    output:
    path("merged_ADO_${label}.tsv"),          emit: merged_ADO
    path("ADO_plot_summary_${label}.png"),    emit: plot_ADO
    path("merged_ADO_summary_${label}.tsv"),  emit: summary_ADO

    script:
    """
    set -euo pipefail

    echo -e "File_Interval\\tFreq\\tProp" > merged_ADO_${label}.tsv

    cat res_ADO_* | cut -f1  | sed 's/\\.tsv//' > a.txt
    cat res_ADO_* | sed 's/\\[//' | sed 's/)//' | sed 's/]//' | cut -f2 > b.txt
    cat res_ADO_* | cut -f3 > c.txt
    cat res_ADO_* | cut -f4 > d.txt

    paste -d "_" a.txt b.txt > end.txt
    paste -d "\\t" end.txt c.txt d.txt >> merged_ADO_${label}.tsv

    Rscript /usr/local/bin/plot_summary_ado_intervals.R merged_ADO_${label}.tsv

    # R script writes merged_ADO_summary.tsv and ADO_plot_summary.png with hardcoded names
    mv merged_ADO_summary.tsv merged_ADO_summary_${label}.tsv
    mv ADO_plot_summary.png   ADO_plot_summary_${label}.png
    """
}
