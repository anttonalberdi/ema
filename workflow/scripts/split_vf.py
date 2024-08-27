import re
import sys
from collections import defaultdict

def split_fasta_by_vf(input_file):
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

    # Write the content to separate files
    for vf_pattern, lines in vf_content.items():
        output_file = f"{vf_pattern}.fasta"
        with open(output_file, 'w') as outfile:
            outfile.writelines(lines)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python workflow/scripts/split_vf.py input_fasta_file")
        sys.exit(1)

    input_file = sys.argv[1]
    split_fasta_by_vf(input_file)
