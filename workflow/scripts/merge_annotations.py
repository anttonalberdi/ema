import pandas as pd
import argparse

def process_files(gff_file, pfam_file, output_file):
    # Read the GFF and PFAM files
    gff_df = pd.read_csv(gff_file, sep='\t', header=None)
    pfam_df = pd.read_csv(pfam_file, sep='\t', header=None)
    
    # Extract the relevant columns (2nd and 3rd columns, which are 0-indexed as 1 and 2)
    gff_cols = gff_df.iloc[:, [1, 2]]
    pfam_cols = pfam_df.iloc[:, [1, 2]]
    
    # Create a new column in GFF DataFrame to store the comparison results
    gff_df['Comparison_Result'] = gff_cols.apply(tuple, axis=1).isin(pfam_cols.apply(tuple, axis=1))
    
    # Replace True/False with the desired values
    gff_df['Comparison_Result'] = gff_df['Comparison_Result'].replace({True: 'Match', False: 'No Match'})
    
    # Write the modified GFF to a new file
    gff_df.to_csv(output_file, sep='\t', index=False, header=False)

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Compare GFF and PFAM files and output the results.')
    parser.add_argument('gff_file', type=str, help='Path to the GFF file')
    parser.add_argument('pfam_file', type=str, help='Path to the PFAM file')
    parser.add_argument('output_file', type=str, help='Path to the output file')
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Process the files
    process_files(args.gff_file, args.pfam_file, args.output_file)

if __name__ == '__main__':
    main()
