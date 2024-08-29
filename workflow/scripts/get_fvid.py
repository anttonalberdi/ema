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
                # Extract the whole line and the first word
                first_word = line.split()[0]
                
                # Use regex to find matches for VF or VFC patterns
                matches = re.findall(r'^>\S+|\bVF\d{4}\b|\bVFC\d{4}\b', line)
                if matches:
                    # Store the first word and matches
                    records.append([first_word.lstrip('>')] + matches)

    # Convert records to a DataFrame for easier handling
    df = pd.DataFrame(records)

    # Write the DataFrame to the output file, without an index and tab-separated
    df.to_csv(output_file, sep='\t', header=False, index=False)

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
