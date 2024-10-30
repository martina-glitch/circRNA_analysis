import csv
import os

# Outputs a CSV file for each circRNA identifier with the amount of expression,
# and sorts the entries in descending order.

# Ask the user to input the path of the input file
filename = input("Enter the input file path: ")

# Extract the base file name (without the path)
base_filename = os.path.basename(filename)

# Initialize a list to store the data to write to the CSV file
data = []

# Variable for counting the total number of circRNAs
total_circrna = 0

# Variable for counting the total number of BSJ (backspliced junctions)
total_bsj = 0

# Open the output file in write mode
output_filename = f'total_bsj_counts_illumina_{base_filename}.csv'
with open(output_filename, 'w', newline='') as csvfile:
    # Create a CSV writer object
    csv_writer = csv.writer(csvfile)
    
    # Write the CSV header
    csv_writer.writerow(['circRNA ID', 'BSJ'])

    # Open the input file and read it line by line
    with open(filename, 'r') as file:
        # Skip the header line
        next(file)
        
        # Iterate over each line in the file
        for line in file:
            # Split the line into fields using tab as the delimiter
            fields = line.strip().split('\t')
            
            # Extract the circRNA ID and BSJ count from the line
            circrna_id = fields[0].strip()
            junctions = int(fields[4])
            
            # Add the circRNA ID and BSJ count to the list
            data.append((circrna_id, junctions))
            
            # Update the total counts
            total_circrna += 1
            total_bsj += junctions

    # Sort the list in descending order based on the BSJ count
    data.sort(key=lambda x: x[1], reverse=True)

    # Write the sorted data to the CSV file
    for circrna_id, junctions in data:
        csv_writer.writerow([circrna_id, junctions])

# Print the total number of circRNAs and total BSJ count
print(f"Total number of circRNAs: {total_circrna}")
print(f"Total number of BSJs: {total_bsj}")

# Print a confirmation message with the correct output file name
print(f"Output saved to {output_filename}")


# old code
"""
# For each different Illumina file (A1, A2, ..., A6)
# It takes output data from CIRI2 (e.g., outfile_ciri2_fastp_A1)
# It distinguishes by identifier and counts how much each circRNA is expressed

# Ask the user to input the path of the input file
filename = input("Enter the input file path: ")

# Initialize a counter for backspliced junctions
junction_count = 0
# Initialize a set to store unique circRNA IDs
circrna_ids = set()

# Open the file and read it line by line
with open(filename, 'r') as file:
    # Skip the header line
    next(file)
    # Iterate over each line in the file
    for line in file:
        # Split the line into fields using tab as the delimiter
        fields = line.split('\t')
        # Extract the number of backspliced junctions from the appropriate column
        # In this case, the fifth column (index 4)
        junctions = int(fields[4])
        # Increment the total count of backspliced junctions
        junction_count += junctions

        circrna_id = fields[0].strip()
        circrna_ids.add(circrna_id)

# Get the total number of different types of circRNA by counting unique IDs
total_circrna_types = len(circrna_ids)

# Print the total number of different types of circRNA
print("Total number of different types of circRNA:", total_circrna_types)

# Print the total number of backspliced junctions, representing the number of circRNAs
print("Total number of circRNAs (backspliced junctions):", junction_count)
"""
