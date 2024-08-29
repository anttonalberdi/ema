import re
import pandas as pd
import sys

def get_fvid(input_file, output_file):
    # Initialize a list to store extracted data
    records = []

    # Open the input file and process each line
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Extract the first word (e.g., VFG037170(gb|WP_001081754))
                first_word = line.split()[0].lstrip('>')
                
                # Use regex to find matches for VF or VFC patterns
                vf_matches = re.findall(r'\bVF\d{4}\b', line)
                vfc_matches = re.findall(r'\bVFC\d{4}\b', line)
                
                # Extract the text between ") - " and "(VFC"
                match = re.search(r'\) - (.*?)(?:\(VFC|$)', line)
                description = match.group(1).strip() if match else ''
                
                # Append the extracted data to the records list
                records.append([first_word] + vf_matches + vfc_matches + [description])

    # Convert records to a DataFrame for easier handling
    df = pd.DataFrame(records, columns=['entry', 'vf', 'vfc', 'vf_type'])

    # Write the DataFrame to the output file, without an index and tab-separated
    df.to_csv(output_file, sep='\t', header=True, index=False)

if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python script_name.py <input_file> <output_file>")
        sys.exit(1)

    # Assign command-line arguments to variables
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Call the processing function
    get_fvid(input_file, output_file)
