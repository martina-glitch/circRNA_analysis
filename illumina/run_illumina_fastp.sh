#!/bin/bash

# Script for filtering Illumina samples using fastp
# This script processes Illumina samples A2, A3, and A4, performing adapter trimming,
# error correction, and quality filtering using fastp. The output includes both filtered
# sequence files and an HTML quality report for each sample.

# List of samples to process
samples=("A2" "A3" "A4")

# Define input and output directories
input_path="/data/mariachiara/illumina_circRNA/X204SC23075032-Z01-F001_02/01.RawData"
output_path="/data/martina.zambon/results_illumina_fastp"

# Iterate over each sample in the list
for index in "${!samples[@]}"; do
    sample="${samples[$index]}"

    # Construct input file paths for paired-end reads
    input_forward="${input_path}/${sample}/${sample}_FKRN24001880$((index))-1A_H7C75DSXC_L4_1.fq.gz"
    input_reverse="${input_path}/${sample}/${sample}_FKRN24001880$((index))-1A_H7C75DSXC_L4_2.fq.gz"

    # Define output file prefix and HTML report path
    output_prefix="${output_path}/out_${sample}"
    html_output="${output_path}/fastp_output_${sample}.html"
    
    # Build the fastp command
    fastp_command=(
        "fastp"
        "-i" "$input_forward"                   # Input file for read 1
        "-o" "${output_prefix}_1"               # Output file for read 1 after filtering
        "-I" "$input_reverse"                   # Input file for read 2
        "-O" "${output_prefix}_2"               # Output file for read 2 after filtering
        "${output_prefix}_unpaired1"            # Output file for unpaired reads from read 1
        "${output_prefix}_unpaired2"            # Output file for unpaired reads from read 2
        "${output_prefix}_failed_out"           # Output file for reads that fail quality filtering
        "${output_prefix}_overlapped_out"       # Output file for overlapped reads
        "-c" "--correction"                     # Enable error correction
        "--overrepresentation_analysis"         # Analyze overrepresented sequences
        "--adapter_sequence" "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"    # Adapter for read 1
        "--adapter_sequence_r2" "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG" # Adapter for read 2
        "-D" "--dedup"                          # Perform duplicate removal
        "--dup_calc_accuracy" "6"               # Set duplicate calculation accuracy level
        "--n_base_limit" "3"                    # Limit the number of ambiguous N bases
        "-h" "$html_output"                     # Generate an HTML report
        "--report_title" "fast_report_${sample}" # Title for the report
        "--thread" "8"                          # Number of threads to use
        "--trim_poly_g"                         # Trim poly-G tails
        "--trim_poly_x"                         # Trim poly-X tails
        "--low_complexity_filter"               # Filter out low-complexity sequences
    )

    # Execute the fastp command
    "${fastp_command[@]}"
done
