import subprocess
import os
import glob
import logging

# This script processes Nanopore sequencing data from .fq.gz files by performing the following steps:
# 1. Filters the sequencing reads using NanoFilt with a minimum quality score and length threshold.
# 2. Converts the filtered reads from FASTQ format to FASTA format.
# 3. Maps the filtered reads to a reference genome using pblat, generating PSL alignment files.
# 4. Logs any errors encountered during filtering or mapping.
# 5. Concatenates all filtered FASTA files into a single output file.
# Optionally, it can also convert PSL files to SAM format for further analysis (this part is currently commented out).

def filter_and_map(directory, fa):
    total_reads = 0
    fasta_files = []
    psl_files = []
    output_dir = "/data/martina.zambon/results_nanopore_filtered_9"
    os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists
    
    # Find all .fq.gz files in the specified directory
    file_list = glob.glob(os.path.join(directory, "*.fq.gz"))

    for file_name in file_list:
        output_fasta = os.path.join(output_dir, f"{os.path.basename(file_name)}_filtered.fa")
        
        # Use NanoFilt to filter reads and convert them to FASTA format
        command = (
            f"zcat {file_name} | NanoFilt -q 7 -l 200 | "
            f"sed -n '1~4s/^@/>/p;2~4p' > {output_fasta}"
        )

        # Execute the command in the shell
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        
        if process.returncode == 0:
            # Count the number of reads (each read has 2 lines in FASTA)
            total_reads += len(open(output_fasta).readlines()) // 2
            fasta_files.append(output_fasta)  # Save the filtered FASTA file
        else:
            # Log any errors encountered during NanoFilt execution
            logging.error(f"Error executing NanoFilt for file {file_name}: {stderr.decode().strip()}")

        # Run pblat to align the filtered reads and save the output to a PSL file
        psl_file = os.path.join(output_dir, f"{os.path.basename(file_name)}.psl")
        psl_files.append(psl_file)
        pblat_command = f"pblat -threads=8 -trimT -mask=lower {fa} {output_fasta} {psl_file}"
        
        # Execute the pblat command and capture any errors
        pblat_process = subprocess.run(pblat_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if pblat_process.returncode != 0:
            logging.error(f"Error executing PBLAT for file {file_name}: {pblat_process.stderr.decode().strip()}")

    # If needed, convert PSL files to SAM format (this part is currently commented out)
    # sam_output_file = os.path.join(output_dir, "output.sam")
    # with open(sam_output_file, "w") as sam_file:
    #     for psl_file in psl_files:
    #         convert_psl_to_sam(psl_file, sam_file)
    
    # Concatenate all the filtered FASTA files into one output file
    with open(os.path.join(output_dir, "output.fasta"), "w") as output_file:
        for fasta_file in fasta_files:
            with open(fasta_file, "r") as input_file:
                output_file.write(input_file.read())
    
    logging.info("Total number of reads after NanoFilt filtering: %d", total_reads)

# Optional function to convert PSL files to SAM format
def convert_psl_to_sam(psl_file, sam_file):
    command = f"uncle_psl.py {psl_file}"
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    sam_output = p.communicate()[0].decode()
    with open(sam_file, "a") as sam_output_file:
        sam_output_file.write(sam_output)

# Set up logging to track the process
logging.basicConfig(level=logging.INFO)

# Path to the directory containing .fq.gz files and reference genome
directory = "/data/martina.zambon/results_nanopore_9_from_fq.gz/"
fa = "/data/martina.zambon/data/human/hg19.fa"

# Run the filtering and mapping process
filter_and_map(directory, fa)
