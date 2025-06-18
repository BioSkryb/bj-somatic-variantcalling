nextflow.enable.dsl=2
params.timestamp = ""

process CREATE_TAB_NVNR {
    tag "${group}_${chr}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"
    
    input:
    tuple val(group),val(chr), path(df_nv), path(df_nr_files)
    val(publish_dir)
    val(enable_publish)
  
    output:
    tuple val(group), val(chr), path("mat_nv_group_${group}_chr_${chr}.tsv"), path("mat_nr_group_${group}_chr_${chr}.tsv")

    script:
    """

    echo -e "Listing content to process ...";

    ls;

    echo -e "Cleaning NV ...";

    head -n1 ${df_nv} > header_nv.tsv

    cat ${df_nv} | grep "^${chr}_" | tail -n+2 | awk -v OFS="\\t" '{outstring=\$1;for(i=2; i<=NF; i=i+1){gsub (/.*,/,"",\$i);gsub(/\\./,"0",\$i);outstring=outstring"\\t"\$i;}print outstring;}'  > df_nv_clean.tsv

    cat header_nv.tsv df_nv_clean.tsv >  mat_nv_group_${group}_chr_${chr}.tsv

    echo -e "Looping of NR files and placing in same order as NV ...";

    cat mat_nv_group_${group}_chr_${chr}.tsv | tail -n+2 | cut -f1 > order_variants_ids_nv.txt

    cat order_variants_ids_nv.txt | cut -d "_" -f1,2 > order_pos_nv.txt

    cat mat_nv_group_${group}_chr_${chr}.tsv| head -n1 | tr "\\t" "\\n" | tail -n+2 > order_samples_nv.txt

    cp order_variants_ids_nv.txt df_nr_clean.tsv

    cat order_samples_nv.txt | while read nom;
    do

        touch temp_nr.tsv;

        ls \${nom}_df_nr_group*.tsv  | sort -V  | while read nr; do cat \${nr} >> temp_nr.tsv;done

        echo -e "Pocessing \${nom}";

        awk -v OFS="\\t" 'NR == FNR {  a[\$1]=\$2; next }{if( \$0 in a ){print a[\$0]}else{print 0}}' temp_nr.tsv order_pos_nv.txt > temp.txt

        paste -d "\\t" df_nr_clean.tsv temp.txt > temp_paste.txt

        mv temp_paste.txt df_nr_clean.tsv

        rm temp.txt

        rm temp_nr.tsv;
    
    done

    cat header_nv.tsv df_nr_clean.tsv >  mat_nr_group_${group}_chr_${chr}.tsv

    wc -l mat*.tsv

    """
}