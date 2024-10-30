import os

# This script processes circRNA output files from the Nanopore tool, merges the results into a single file,
# groups circRNAs by the same ID (same gene) and BSJ (backspliced junction) count,
# calculates how many times each circRNA is repeated, and sorts them in descending order.

# List of directories and files
CIRCRNA_ANNOTATION_DIRS = [
    "/data/martina.zambon/results_nanopore_1/FAY73732_pass_7f61b213_4811f00b_all/", # A1
    "/data/martina.zambon/results_nanopore_2/FAY49674_pass_cedafbe6_cdba8e78_all/", # A2
    "/data/martina.zambon/results_nanopore_3/FAY72835_pass_fc30eac1_a43abb16_all/", # A3
    "/data/martina.zambon/results_nanopore_4/FAY48470_pass_c8aa96a2_7f98c47e_all/", # A4
    "/data/martina.zambon/results_nanopore_5/FAY75373_pass_all/", # A5
    "/data/martina.zambon/results_nanopore_6/FAY74971_pass_92159db0_ec7149eb_all/", # A6
    "/data/martina.zambon/results_nanopore_7/all_A7/",
    "/data/martina.zambon/results_nanopore_8/all_A8/",
    "/data/martina.zambon/results_nanopore_9/all_A9/"
]

OUTPUT_FILES = [
    "outfile_nanopore_control_A1.txt",
    "outfile_nanopore_control_A2.txt",
    "outfile_nanopore_control_A3.txt",
    "outfile_nanopore_control_A4.txt",
    "outfile_nanopore_control_A5.txt",
    "outfile_nanopore_control_A6.txt",
    "outfile_nanopore_control_A7.txt",
    "outfile_nanopore_control_A8.txt",
    "outfile_nanopore_control_A9.txt"
]

CIRC_RNA_READS_FILES = [
    "circRNA_reads_A1.txt",
    "circRNA_reads_A2.txt",
    "circRNA_reads_A3.txt",
    "circRNA_reads_A4.txt",
    "circRNA_reads_A5.txt",
    "circRNA_reads_A6.txt",
    "circRNA_reads_A7.txt",
    "circRNA_reads_A8.txt",
    "circRNA_reads_A9.txt"
]

# Function to process a single sample
def process_sample(annotation_dir, output_file, reads_file):
    # Dictionary to track the total count of BSJ for each RNA identifier
    bsj_counts = {}

    # Analyze the input files
    for dirpath, dirnames, filenames in os.walk(annotation_dir):
        for file in filenames:
            if file.endswith(".circRNA_candidates.annotated.txt"):
                input_file = os.path.join(dirpath, file)

                with open(input_file, 'r') as f_in:
                    # Skip the first row (file header)
                    next(f_in)

                    for line in f_in:
                        parts = line.strip().split('\t')
                        rna_id = parts[0].split("_")[3]  # circRNA ID
                        chr = parts[1]     # Chromosome
                        start = parts[2]   # Start position
                        end = parts[3]     # End position
                        bsj_reads = int(parts[5])  # BSJ count
                        key = (rna_id, chr, start, end)  # Dictionary key
                        bsj_counts[key] = bsj_counts.get(key, 0) + bsj_reads

    # Calculate the total number of BSJ
    total_bsj = sum(bsj_counts.values())

    # Calculate the total number of distinct circRNA IDs
    num_distinct_rna_id = len(set([key for key in bsj_counts.keys()]))

    # Print the total number of BSJ and distinct circRNA IDs
    print(f"Sample: {annotation_dir}")
    print("Total BSJ count:", total_bsj)
    print("Total distinct circRNA IDs:", num_distinct_rna_id)

    # Sort the dictionary by BSJ count in descending order
    sorted_bsj_counts = sorted(bsj_counts.items(), key=lambda x: x[1], reverse=True)

    # Read the input data from the circRNA_reads.txt file
    input_data = []
    with open(reads_file, 'r') as f_input:
        next(f_input)  # Skip the header
        for line in f_input:
            parts = line.strip().split('\t')
            rna_id = parts[0]
            chr = parts[1]
            start = parts[2]
            end = parts[3]
            id_read = parts[4].split('~')[0]  # Take only the first part of ID_read before '~'
            circBase = parts[6]
            circAtlas = parts[7]
            circpedia = parts[8]
            input_data.append((rna_id, chr, start, end, id_read, circBase, circAtlas, circpedia))

    # Write the output data to the output file
    with open(output_file, 'w') as f_out:
        f_out.write("circRNA_name\tchr\tstart\tend\tID_read\tbsj_reads\tcircBase_ID\tcircAtlas_ID\tcircpedia_ID\tnumero_bsj\n")
        for key, bsj_count in sorted_bsj_counts:
            rna_id, chr, start, end = key
            id_read = ""  # Initialize ID_read as an empty string
            circBase = ""
            circAtlas = ""
            circpedia = ""
            # Find the corresponding ID_read for the key-tuple in the input data
            for rna_id_input, chr_input, start_input, end_input, id_read_input, circBase_input, circAtlas_input, circpedia_input in input_data:
                if rna_id_input == rna_id and chr_input == chr and start_input == start and end_input == end:
                    id_read = id_read_input  # Assign the corresponding ID_read if found
                    circBase = circBase_input
                    circAtlas = circAtlas_input
                    circpedia = circpedia_input
                    break
            f_out.write(f"{rna_id}\t{chr}\t{start}\t{end}\t{id_read}\t{bsj_count}\t{circBase}\t{circAtlas}\t{circpedia}\t{bsj_count}\n")

    print(f"Output completed for sample: {output_file}")

# Process all samples
for i, annotation_dir in enumerate(CIRCRNA_ANNOTATION_DIRS):
    process_sample(annotation_dir, OUTPUT_FILES[i], CIRC_RNA_READS_FILES[i])
