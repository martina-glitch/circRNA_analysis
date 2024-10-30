#!/bin/bash

# Definisce il percorso della cartella contenente i file
# input_folder="/data/mariachiara/circRNA2_nanopore/"
#input_folder="/data/martina.zambon/results_nanopore_1_from_fq.gz/"

input_folder="/data/martina.zambon/results_nanopore_1_from_fq.gz/"


filename="${input_folder}FAY73732_pass_7f61b213_4811f00b_all.fq.gz"
# filename="${input_folder}all_A7.fq.gz"
echo "Eseguo long_read_circRNA su ${filename}..."
long_read_circRNA run "${filename}" -o /data/martina.zambon/results_nanopore_1

# long_read_circRNA run /data/martina.zambon/results_nanopore_1_from_fq.gz/FAY73732_pass_7f61b213_4811f00b_all.fq.gz -o /data/martina.zambon/results_nanopore_1
