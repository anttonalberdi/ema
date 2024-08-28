import argparse
import pandas as pd
from collections import defaultdict
from Bio import SearchIO

def process_files(gff_file, kofams_file, pfam_file, cazy_file, output_file):

    ##############
    # Load genes #
    ##############

    annotations = pd.read_csv(gff_file, sep='\t', comment='#', header=None, 
                     names=['seqid', 'source', 'type', 'start', 'end', 
                            'score', 'strand', 'phase', 'attributes'])
    annotations = annotations.drop(columns=['attributes', 'source', 'score', 'type', 'phase'])
    annotations = annotations.rename(columns={'seqid': 'gene'})

    #####################
    # Parse annotations #
    #####################
    
    hmm_attribs = ['bitscore', 'evalue', 'id', 'overlap_num', 'region_num']

    # Parse KOFAMS
    kofams_hits = defaultdict(list)
    query_ids = []
    
    with open(kofams_file) as handle:
        for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
            query_id = queryresult.id  # Capture the query result id
            query_ids.extend([query_id] * len(queryresult.hits))  # Extend query_ids list to match the number of hits
            
            for hit in queryresult.hits:
                for attrib in hmm_attribs:
                    # Use `getattr` to fetch the attribute value from the hit object
                    value = getattr(hit, attrib, None)  # Use `None` as a default if the attribute does not exist
                    kofams_hits[attrib].append(value)
    
    kofams_hits['query_id'] = query_ids
    kofams_df = pd.DataFrame.from_dict(kofams_hits)
    kofams_df = kofams_df.rename(columns={'query_id': 'gene'})
    
    # Parse PFAM
    pfam_hits = defaultdict(list)
    query_ids = []
    
    with open(pfam_file) as handle:
        for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
            query_id = queryresult.id  # Capture the query result id
            query_ids.extend([query_id] * len(queryresult.hits))  # Extend query_ids list to match the number of hits
            
            for hit in queryresult.hits:
                for attrib in hmm_attribs:
                    # Use `getattr` to fetch the attribute value from the hit object
                    value = getattr(hit, attrib, None)  # Use `None` as a default if the attribute does not exist
                    pfam_hits[attrib].append(value)
    
    pfam_hits['query_id'] = query_ids
    pfam_df = pd.DataFrame.from_dict(pfam_hits)
    pfam_df = pfam_df.rename(columns={'query_id': 'gene'})
    
    # Parse CAZY
    cazy_hits = defaultdict(list)
    query_ids = []
    
    with open(cazy_file) as handle:
        for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
            query_id = queryresult.id  # Capture the query result id
            query_ids.extend([query_id] * len(queryresult.hits))  # Extend query_ids list to match the number of hits
            
            for hit in queryresult.hits:
                for attrib in hmm_attribs:
                    # Use `getattr` to fetch the attribute value from the hit object
                    value = getattr(hit, attrib, None)  # Use `None` as a default if the attribute does not exist
                    cazy_hits[attrib].append(value)
    
    cazy_hits['query_id'] = query_ids
    cazy_df = pd.DataFrame.from_dict(cazy_hits)
    cazy_df = cazy_df.rename(columns={'query_id': 'gene'})

    # Parse AMR
    amr_hits = defaultdict(list)
    query_ids = []
    
    with open(amr_file) as handle:
        for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
            query_id = queryresult.id  # Capture the query result id
            query_ids.extend([query_id] * len(queryresult.hits))  # Extend query_ids list to match the number of hits
            
            for hit in queryresult.hits:
                for attrib in hmm_attribs:
                    # Use `getattr` to fetch the attribute value from the hit object
                    value = getattr(hit, attrib, None)  # Use `None` as a default if the attribute does not exist
                    amr_hits[attrib].append(value)
    
    amr_hits['query_id'] = query_ids
    amr_df = pd.DataFrame.from_dict(amr_hits)
    amr_df = amr_df.rename(columns={'query_id': 'gene'})

    
    
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
    parser.add_argument('-gff', required=True, type=str, help='Path to the GFF file')
    parser.add_argument('-pfam', required=True, type=str, help='Path to the PFAM file')
    parser.add_argument('-o', required=True, type=str, help='Path to the output file')
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Process the files
    process_files(args.gff, args.pfam, args.o)

if __name__ == '__main__':
    main()
