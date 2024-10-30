# Normalizes the expression of circRNAs from Illumina by the number of reads with BSJ
# and multiplies it by one million.

import csv

# Function to normalize expression values
def normalize_expression(values, total_reads):
    normalized_values = []
    for value, total in zip(values, total_reads):
        if total == 0:
            normalized_values.append(0.0)
        else:
            normalized_value = (float(value) / total) * 1000000
            normalized_values.append(normalized_value)
    return normalized_values

# Reading the input CSV file
input_filename = 'sorted_Amerge1.csv'
output_filename = 'normalized_sorted_Amerge1.csv'

# Step 1: Calculate the total number of reads for each sample
total_reads_per_sample = None

with open(input_filename, 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=';')
    
    for row in reader:
        if len(row) < 4:
            continue  # If the row has fewer than 4 elements, skip it
        
        samples = list(map(int, row[3:]))
        
        if total_reads_per_sample is None:
            total_reads_per_sample = samples
        else:
            total_reads_per_sample = [total + sample for total, sample in zip(total_reads_per_sample, samples)]

# Step 2: Normalize the expression values and write to the output file
with open(input_filename, 'r') as csvfile, open(output_filename, 'w', newline='') as outfile:
    reader = csv.reader(csvfile, delimiter=';')
    writer = csv.writer(outfile, delimiter=';')
    
    for row in reader:
        if len(row) < 4:
            continue  # If the row has fewer than 4 elements, skip it
        
        chromosome = row[0]
        start = row[1]
        end = row[2]
        samples = list(map(int, row[3:]))
        
        # Normalize the expression values
        normalized_samples = normalize_expression(samples, total_reads_per_sample)
        
        # Build the output row
        output_row = [chromosome, start, end] + normalized_samples
        
        # Write the row to the output file
        writer.writerow(output_row)

print(f"File '{output_filename}' successfully created.")
