import os

# Path of the input file
input_file_path = '/data/martina.zambon/results_illumina_fastp/out_A1_1'

# Path of the output file
output_file_path = 'count_reads_illumina.txt'

# Check if the input file path exists
if not os.path.exists(input_file_path):
    print("The input file does not exist.")
    exit()

# Open the input file
with open(input_file_path, 'r') as input_file:
    lines = input_file.readlines()

# Open the output file
with open(output_file_path, 'w') as output_file:
    # Iterate over every block of 4 lines (identifier, sequence, "+", and quality)
    for i in range(0, len(lines), 4):
        identifier = lines[i].strip()  # First line: identifier
        sequence = lines[i + 1].strip()  # Second line: sequence
        length = len(sequence)  # Calculate the length of the sequence
        
        # Write the length to the output file
        output_file.write(f'{length}\n')

print("Read length counting completed. The result has been saved in", output_file_path)
