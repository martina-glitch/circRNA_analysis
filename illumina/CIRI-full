#!/bin/bash

# This script processes a list of Illumina samples (A1, A2, ..., A9) through multiple steps:
# 1. Indexes the reference genome using BWA.
# 2. Trims reads to a uniform length for each sample.
# 3. Repairs the trimmed paired-end reads to ensure correct pairing.
# 4. Aligns the reads using BWA MEM.
# 5. Runs CIRI2 to detect circRNAs.
# 6. Runs CIRI-AS for circRNA splicing analysis.
# 7. Runs CIRI-full for full-length circRNA reconstruction.
# 8. Visualizes the results using CIRI-vis.

# Define the list of samples to process
samples=("A1" "A2" "A3" "A4" "A5" "A6" "A7" "A8" "A9")

# Define paths for input and output files
output_path="/data/martina.zambon/results_illumina_fastp"
ref_genome="/data/martina.zambon/data/human/hg19.fa"
ref_gtf="/data/martina.zambon/data/human/hg19.refGene.gtf"
trimmed_output_suffix="_trimmed.fastq"
length=150  # Define the uniform trimming length

# Paths to scripts
ciri2_script="/home/martina.zambon/CIRI_v2.0.6/CIRI2.pl"
ciri_as_script="/home/martina.zambon/CIRI_AS.pl"
ciri_full_jar="/home/martina.zambon/CIRI-full_v2.0/CIRI-full.jar"
ciri_vis_jar="/home/martina.zambon/CIRI-vis.jar"
bbmap_dir="/home/martina.zambon/bbmap"

# Index the reference genome using BWA
echo "Indexing the reference genome..."
bwa index -a bwtsw $ref_genome

# Trim the reads to a uniform length
for sample in "${samples[@]}"; do
    input_r1="${output_path}/out_${sample}_1"
    input_r2="${output_path}/out_${sample}_2"
    trimmed_r1="${output_path}/trimmed_out_${sample}_1${trimmed_output_suffix}"
    trimmed_r2="${output_path}/trimmed_out_${sample}_2${trimmed_output_suffix}"

    echo "Trimming reads for sample ${sample}..."
    cutadapt -j 10 -l $length -m $length \
             -o "${trimmed_r1}" -p "${trimmed_r2}" \
             "${input_r1}" "${input_r2}"
    
    # Repair trimmed reads to ensure paired-end reads are correctly matched
    echo "Repairing reads for sample ${sample}..."
    $bbmap_dir/repair.sh in1="${trimmed_r1}" in2="${trimmed_r2}" \
                         out1="${trimmed_r1}" out2="${trimmed_r2}" \
                         outs="${output_path}/singletons_${sample}.fastq"
done

# Align the trimmed and repaired reads using BWA MEM
for sample in "${samples[@]}"; do
    trimmed_r1="${output_path}/trimmed_out_${sample}_1${trimmed_output_suffix}"
    trimmed_r2="${output_path}/trimmed_out_${sample}_2${trimmed_output_suffix}"
    alignment_output="${output_path}/alignment_${sample}_aligned.sam"

    echo "Aligning reads for sample ${sample}..."
    bwa mem -t 8 -T 19 $ref_genome "${trimmed_r1}" "${trimmed_r2}" > "${alignment_output}"
done

# Run CIRI2 for circRNA detection
for sample in "${samples[@]}"; do
    alignment_output="${output_path}/alignment_${sample}_aligned.sam"
    ciri_output="${output_path}/outfile_ciri2_${sample}_fastp"

    echo "Running CIRI2 for sample ${sample}..."
    perl $ciri2_script \
        -I "${alignment_output}" \
        -O "${ciri_output}" \
        -F $ref_genome \
        -A $ref_gtf \
        -T 10
done

# Run CIRI-AS for splicing analysis
for sample in "${samples[@]}"; do
    alignment_output="${output_path}/alignment_${sample}_aligned.sam"
    ciri_output="${output_path}/outfile_ciri2_${sample}_fastp"
    ciriAS_output="${output_path}/outfile_ciriAS_${sample}_fastp"

    # Check if the CIRI-AS output file exists and delete it if necessary
    if [ -f "${ciriAS_output}.list" ]; then
        echo "File ${ciriAS_output}.list already exists, deleting it."
        rm "${ciriAS_output}.list"
    fi

    echo "Running CIRI-AS for sample ${sample}..."
    perl $ciri_as_script \
        -S "${alignment_output}" \
        -C "${ciri_output}" \
        -O "${ciriAS_output}" \
        -F $ref_genome \
        -A $ref_gtf \
        -D yes
done

# Run CIRI-full for full-length circRNA reconstruction for each sample
for sample in "${samples[@]}"; do
    trimmed_r1="${output_path}/trimmed_out_${sample}_1${trimmed_output_suffix}"
    trimmed_r2="${output_path}/trimmed_out_${sample}_2${trimmed_output_suffix}"
    prefix="${output_path}/prefix_${sample}"
    output_dir="${output_path}/CIRI_full_output_${sample}"

    # CIRI-full Pipeline
    echo "Running CIRI-full Pipeline for sample ${sample}..."
    java -jar $ciri_full_jar Pipeline \
        -1 "${trimmed_r1}" \
        -2 "${trimmed_r2}" \
        -r $ref_genome \
        -a $ref_gtf \
        -o "prefix_${sample}" \
        -d "${output_dir}" \
        -t 4  # Use 4 threads

    # CIRI-full RO1
    echo "Running CIRI-full RO1 for sample ${sample}..."
    java -jar $ciri_full_jar RO1 \
        -1 "${trimmed_r1}" \
        -2 "${trimmed_r2}" \
        -o "prefix_${sample}_ro1"

    # BWA MEM alignment for RO1 output
    bwa_mem_output="${prefix}_ro1.sam"
    bwa_mem_input="${prefix}_ro1.fq"
    echo "Aligning reads for RO1 of sample ${sample}..."
    bwa mem -T 19 $ref_genome "${bwa_mem_input}" > "${bwa_mem_output}"

    # CIRI-full RO2
    echo "Running CIRI-full RO2 for sample ${sample}..."
    java -jar $ciri_full_jar RO2 \
        -r $ref_genome \
        -s "${bwa_mem_output}" \
        -l $length \
        -o "prefix_${sample}"

    # CIRI-full Merge
    ciri2_output="${output_path}/outfile_ciri2_${sample}_fastp"
    ciriAS_output="${output_path}/outfile_ciriAS_${sample}_fastp"
    ro2_info_list="${prefix}_ro2_info.list"
    echo "Running CIRI-full Merge for sample ${sample}..."
    java -jar $ciri_full_jar Merge \
        -a $ref_gtf \
        -c "${ciri2_output}" \
        -as "${ciriAS_output}_jav.list" \
        -ro "${ro2_info_list}" \
        -o "prefix_${sample}" \
        -r $ref_genome
done

# Visualize results with CIRI-vis
for sample in "${samples[@]}"; do
    anno_file="${output_path}/prefix_${sample}_merge_circRNA_detail.anno"
    lib_length_file="${output_path}/outfile_ciriAS_${sample}_fastp_library_length.list"
    output_vis="${output_path}/CIRI-vis_${sample}"

    echo "Running CIRI-vis for sample ${sample}..."
    java -jar $ciri_vis_jar -i "${anno_file}" -l "${lib_length_file}" -r $ref_genome -d "${output_vis}"
done

echo "All steps completed."
