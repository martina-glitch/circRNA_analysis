#!/bin/bash

# 1. Align Illumina samples (output from fastp) to the reference genome hg19
# 2. CIRI2 detects circRNAs from the alignment files

# Define the list of samples to process
samples=("A1" "A2" "A3" "A4" "A5" "A6")

# Define the input/output file paths
output_path="/data/martina.zambon/results_illumina_fastp"

# Loop through each sample
for index in "${!samples[@]}"; do
    sample="${samples[$index]}"
    output_prefix="${output_path}/out_${sample}"

    # Step 1: Align each sample using BWA with the hg19 reference genome
    bwa mem /data/martina.zambon/data/human/hg19.fa "${output_prefix}_1" "${output_prefix}_2" > "alignment_fastp_${sample}.sam"

    # Step 2: Run CIRI2 to detect circRNAs using the alignment file
    perl_command=(
        "perl"
        "../../home/martina.zambon/CIRI_v2.0.6/CIRI2.pl"
        "-I" "alignment_fastp_${sample}.sam"
        "-O" "outfile_ciri2_fastp_${sample}.txt"
        "-F" "/data/martina.zambon/data/human/hg19.fa"
        "-A" "/data/martina.zambon/data/human/hg19.refGene.gtf"
    )

    # Execute the CIRI2 command
    "${perl_command[@]}"

done
