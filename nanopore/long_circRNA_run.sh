#!/bin/bash

# Define the path to the folder containing the files
# input_folder="/data/mariachiara/circRNA2_nanopore/"
# input_folder="/data/martina.zambon/results_nanopore_1_from_fq.gz/"

input_folder="/data/martina.zambon/results_nanopore_1_from_fq.gz/"

# Define the filename to process
filename="${input_folder}FAY73732_pass_7f61b213_4811f00b_all.fq.gz"
# Alternative filename: 
# filename="${input_folder}all_A7.fq.gz"

# Execute long_read_circRNA tool on the specified file
echo "Running long_read_circRNA on ${filename}..."
long_read_circRNA run "${filename}" -o /data/martina.zambon/results_nanopore_1

# Example command (for reference):
# long_read_circRNA run /data/martina.zambon/results_nanopore_1_from_fq.gz/FAY73732_pass_7f61b213_4811f00b_all.fq.gz -o /data/martina.zambon/results_nanopore_1
