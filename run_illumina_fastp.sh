#!/bin/bash

# filtra i campioni illumina A2 A3 e A4

# Definisci la lista dei campioni da modificare
campioni=("A2" "A3" "A4")

# Percorso dei file di input e output
percorso_input="/data/mariachiara/illumina_circRNA/X204SC23075032-Z01-F001_02/01.RawData"
percorso_output="/data/martina.zambon/results_illumina_fastp"

# Loop attraverso i campioni
for indice in "${!campioni[@]}"; do
    campione="${campioni[$indice]}"
    # Costruisci i nomi dei file di input e output
    input_forward="${percorso_input}/${campione}/${campione}_FKRN24001880$((indice))-1A_H7C75DSXC_L4_1.fq.gz"
    input_reverse="${percorso_input}/${campione}/${campione}_FKRN24001880$((indice))-1A_H7C75DSXC_L4_2.fq.gz"
    output_prefix="${percorso_output}/out_${campione}"
    html_output="${percorso_output}/uscita_fastp_${campione}.html"
    

    # Costruisci il comando fastp
    comando_fastp=(
        "fastp"
        "-i" "$input_forward"
        "-o" "${output_prefix}_1"
        "-I" "$input_reverse"
        "-O" "${output_prefix}_2"
        "${output_prefix}_unpaired1"
        "${output_prefix}_unpaired2"
        "${output_prefix}_failed_out"
        "${output_prefix}_overlapped_out"
        "-c" "--correction"
        "--overrepresentation_analysis"
        "--adapter_sequence" "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
        "--adapter_sequence_r2" "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG"
        "-D" "--dedup" "--dup_calc_accuracy" "6"
        "--n_base_limit" "3"
        "-h" "$html_output" "--report_title" "fast_report_${campione}"
        "--thread" "8"
        "--trim_poly_g"
        "--trim_poly_x"
        "--low_complexity_filter"
        
    )

    # Esegui il comando Fastp
    "${comando_fastp[@]}"
done

