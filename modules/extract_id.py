def extract_id(result, id_list, evaluation_threshold):
    import re
    import os

    file_path = os.path.expanduser(result) 

    output_path = os.path.expanduser(id_list)  

    with open(file_path, 'r') as file:
        data = file.readlines()

    pattern = re.compile('\d+e-[0-9]{2}', re.IGNORECASE)

    process_data = True
    with open(output_path, 'w') as out_file:
        for line in data:
            if '------ inclusion threshold ------' in line:
                process_data = False

            if process_data:
                match = pattern.search(line)
                if match:
                    e_value = float(match.group())

                    if e_value < float(evaluation_threshold) :
                        parts = line.split()
                        out_file.write(parts[-1] + '\n')