import re
import sys
import os
from collections import defaultdict

def split_fasta_by_vf(input_file, output_dir):
    # Dictionary to store the content for each VF pattern
    vf_content = defaultdict(list)

    with open(input_file, 'r') as infile:
        current_vf = None
        for line in infile:
            if line.startswith('>'):
                # Extract the VF pattern
                match = re.search(r'\(VF\d{4}\)', line)
                if match:
                    vf_pattern = match.group(0).replace('(', '').replace(')', '')  # Remove parentheses
                    current_vf = vf_pattern
            
            # Append the line to the corresponding VF pattern
            if current_vf:
                vf_content[current_vf].append(line)

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Write the content to separate files in the specified output directory
    for vf_pattern, lines in vf_content.items():
        output_file = os.path.join(output_dir, f"{vf_pattern}.fasta")
        with open(output_file, 'w') as outfile:
            outfile.writelines(lines)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python workflow/scripts/split_vf.py input_fasta_file output_directory")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    split_vf(input_file, output_dir)
