import re
import os
import input

file_path = os.path.expanduser('./result.out') 

output_path = os.path.expanduser('./id_list.txt')  

with open(file_path, 'r') as file:
    data = file.readlines()

# Pattern to match evalue anywhere in the string
pattern = re.compile('\d+e-[0-9]{2}', re.IGNORECASE)

process_data = True

# Prepare the output file
with open(output_path, 'w') as out_file:
    for line in data:
        if '------ inclusion threshold ------' in line:
            process_data = False

        if process_data:
            match = pattern.search(line)
            if match:
                e_value = float(match.group())

                if e_value < float(input.evaluation_threshold) :
                    parts = line.split()
                    out_file.write(parts[-1] + '\n')  # Write Sequence gene id to file