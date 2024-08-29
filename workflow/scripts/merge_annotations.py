import argparse
import pandas as pd
from collections import defaultdict
from Bio import SearchIO

#######################
# Auxiliary functions #
#######################

# Function to select the row with the lowest evalue or randomly if there is a tie
def select_lowest_evalue(group):
    # Sort the group by 'evalue' and take the first row of the sorted group
    return group.sort_values(by='evalue').head(1)

# Function to select the row with the lowest evalue or randomly if there is a tie
def select_highest_confidence(group):
    # Sort the group by 'evalue' and take the first row of the sorted group
    return group.sort_values(by='confidence', ascending=False).head(1)

# Function to extract the part after the underscore in ID and append it to seqid
def append_suffix_to_seqid(row):
    # Extract the part after 'ID=' and before the first ';'
    id_part = row['attributes'].split(';')[0].replace('ID=', '')
    # Extract the part after the underscore
    suffix = id_part.split('_')[-1]
    return f"{row['seqid']}_{suffix}"

#################
# Main function #
#################

def merge_annotations(gff_file, kofams_file, pfam_file, cazy_file, ec_file, vfdb_file, vf_file, amr_file, signalp_file, output_file):

    ##############
    # Load genes #
    ##############

    annotations = pd.read_csv(gff_file, sep='\t', comment='#', header=None, 
                     names=['seqid', 'source', 'type', 'start', 'end', 
                            'score', 'strand', 'phase', 'attributes'])
    annotations['seqid'] = annotations.apply(append_suffix_to_seqid, axis=1)
    annotations = annotations.drop(columns=['attributes', 'source', 'score', 'type', 'phase'])
    annotations = annotations.rename(columns={'seqid': 'gene'})

    ######################
    # Load mapping files #
    ######################

    #PFAM to EC
    pfam_to_ec = pd.read_csv(ec_file, sep='\t', comment='#', header=0)
    pfam_to_ec = pfam_to_ec[pfam_to_ec['Type'] == 'GOLD']
    pfam_to_ec = pfam_to_ec.rename(columns={'Confidence-Score': 'confidence'})
    pfam_to_ec = pfam_to_ec.rename(columns={'Pfam-Domain': 'pfam'})
    pfam_to_ec = pfam_to_ec.rename(columns={'EC-Number': 'ec'})
    pfam_to_ec['confidence'] = pd.to_numeric(pfam_to_ec['confidence'], errors='coerce')
    pfam_to_ec = pfam_to_ec.groupby('pfam', group_keys=False).apply(select_highest_confidence)

    #Entry to VF
    entry_to_vf = pd.read_csv(vf_file, sep='\t', comment='#', header=0, names=['entry', 'vf', 'vfc'])
    
    #####################
    # Parse annotations #
    #####################
    
    hmm_attribs = ['accession', 'bitscore', 'evalue', 'id', 'overlap_num', 'region_num']
    evalue_threshold=0.001
    
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
    kofams_df['evalue'] = pd.to_numeric(kofams_df['evalue'], errors='coerce')
    kofams_df = kofams_df[kofams_df['evalue'] < evalue_threshold]
    kofams_df = kofams_df.groupby('gene', group_keys=False).apply(select_lowest_evalue)
    kofams_df = kofams_df.rename(columns={'id': 'kegg'})

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
    pfam_df = pfam_df.rename(columns={'accession': 'pfam'})
    pfam_df['pfam'] = pfam_df['pfam'].str.split('.').str[0]
    pfam_df['evalue'] = pd.to_numeric(pfam_df['evalue'], errors='coerce')
    pfam_df = pfam_df[pfam_df['evalue'] < evalue_threshold]
    pfam_df = pfam_df.groupby('gene', group_keys=False).apply(select_lowest_evalue)
    pfam_df = pd.merge(pfam_df, pfam_to_ec[['pfam', 'ec']], on='pfam', how='left')
    
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
    cazy_df['id'] = cazy_df['id'].str.replace('.hmm', '', regex=False)
    cazy_df['evalue'] = pd.to_numeric(cazy_df['evalue'], errors='coerce')
    cazy_df = cazy_df[cazy_df['evalue'] < evalue_threshold]
    cazy_df = cazy_df.groupby('gene', group_keys=False).apply(select_lowest_evalue)
    cazy_df = cazy_df.rename(columns={'id': 'cazy'})

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
    amr_df['evalue'] = pd.to_numeric(amr_df['evalue'], errors='coerce')
    amr_df = amr_df[amr_df['evalue'] < evalue_threshold]
    amr_df = amr_df.groupby('gene', group_keys=False).apply(select_lowest_evalue)
    amr_df = amr_df.rename(columns={'id': 'amr'})

    # Parse VFDB
    vfdb_df = pd.read_csv(vfdb_file, sep='\t', comment='#', header=None, 
                     names=['gene', 'entry', 'identity', 'length', 'mismatches', 
                            'gaps', 'query_start', 'query_end', 'target_start', 'target_end', 'evalue', 'bitscore'])
    vfdb_df['evalue'] = pd.to_numeric(vfdb_df['evalue'], errors='coerce')
    vfdb_df = vfdb_df[vfdb_df['evalue'] < evalue_threshold]
    vfdb_df = vfdb_df.groupby('gene', group_keys=False).apply(select_lowest_evalue)
    vfdb_df = pd.merge(vfdb_df, entry_to_vf[['entry','vf','vfc']], on='entry', how='left')

    # Parse SIGNALP
    signalp_df = pd.read_csv(signalp_file, sep='\t', comment='#', header=None, names=['gene', 'signalp', 'confidence'])
    signalp_df['confidence'] = pd.to_numeric(signalp_df['confidence'], errors='coerce')
    signalp_df = signalp_df.groupby('gene', group_keys=False).apply(select_highest_confidence)
    
    #####################
    # Merge annotations #
    #####################

    annotations = pd.merge(annotations, kofams_df[['gene', 'kegg']], on='gene', how='left')
    annotations = pd.merge(annotations, pfam_df[['gene', 'pfam', 'ec']], on='gene', how='left')
    annotations = pd.merge(annotations, cazy_df[['gene', 'cazy']], on='gene', how='left')
    annotations = pd.merge(annotations, amr_df[['gene', 'amr']], on='gene', how='left')
    annotations = pd.merge(annotations, vfdb_df[['gene', 'vf', 'vfc']], on='gene', how='left')
    annotations = pd.merge(annotations, signalp_df[['gene', 'signalp']], on='gene', how='left')

    # Output the final DataFrame to the output file
    annotations.to_csv(output_file, sep='\t', index=False)

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Compare GFF and PFAM files and output the results.')
    parser.add_argument('-gff', required=True, type=str, help='Path to the GFF file')
    parser.add_argument('-kofams', required=True, type=str, help='Path to the KOFAMS file')
    parser.add_argument('-pfam', required=True, type=str, help='Path to the PFAM file')
    parser.add_argument('-ec', required=True, type=str, help='Path to the EC file')
    parser.add_argument('-cazy', required=True, type=str, help='Path to the CAZY file')
    parser.add_argument('-vfdb', required=True, type=str, help='Path to the VFDB file')
    parser.add_argument('-vf', required=True, type=str, help='Path to the VF file')
    parser.add_argument('-amr', required=True, type=str, help='Path to the AMR file')
    parser.add_argument('-signalp', required=True, type=str, help='Path to the SIGNALP file')
    parser.add_argument('-o', required=True, type=str, help='Path to the OUTPUT file')
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Process the files
    merge_annotations(args.gff, args.kofams, args.pfam, args.cazy, args.ec, args.vfdb, args.vf, args.amr, args.signalp, args.o)

if __name__ == '__main__':
    main()
