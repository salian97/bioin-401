import pandas as pd
import requests
import io
import time
import os
import shutil

PROTEIN_FILEPATH = r'human_enzyme_gpcr_targets_v1.csv'
UNIPROT_API = r'https://rest.uniprot.org/uniprotkb/accessions'
TEMPLATE_KEY_FILEPATH = r'Experimental_Protein_Ligand\0_Experimantal_Proetin_Ligand.xlsx'
TEMPLATES_DIR = r'Experimental_Protein_Ligand'
PROTEIN_OUTPUT_FILEPATH = r'output_enzyme_gpcr_targets_with_sequences.csv'

def fetch_uniprot_sequences(df: pd.DataFrame, id_column='UniProt ID', chunk_size=500):
    """
    Fetches protein sequences from UniProt in batches and maps them to the DataFrame.
    """
    
    # 1. Get unique IDs to avoid redundant API calls
    unique_ids = df[id_column].unique().tolist()
    sequence_map = {}
    
    print(f"Starting sequence retrieval for {len(unique_ids)} unique UniProt IDs...")

    # 2. Process in chunks (UniProt prefers batches over individual hits)
    for i in range(0, len(unique_ids), chunk_size):
        chunk = unique_ids[i:i + chunk_size]
        ids_query = ",".join(chunk)
        
        # parameters for mapping/retrieval
        params = {
            "format": "tsv",
            "accessions": ids_query,
            "fields": "accession,sequence"
        }
        
        try:
            response = requests.get(UNIPROT_API, params=params)
            response.raise_for_status()
            
            # 3. Parse the TSV response into a temporary dict
            # We use sep='\t' because UniProt returns Tab-Separated Values
            batch_df = pd.read_csv(io.StringIO(response.text), sep='\t')
            
            for _, row in batch_df.iterrows():
                sequence_map[row['Entry']] = row['Sequence']
                
            print(f"Fetched {min(i + chunk_size, len(unique_ids))}/{len(unique_ids)}...")
            
            # Be a good citizen: brief pause between batches
            time.sleep(0.5) 
            
        except Exception as e:
            print(f"Error fetching chunk starting at index {i}: {e}")

    # 4. Map the dictionary back to the original DataFrame
    df['Sequence'] = df[id_column].map(sequence_map)
    
    print(f'Fetched {len(unique_ids)} unique protein sequences')
    return df

def mapPDBtemplates(main_df: pd.DataFrame, excel_path, source_pdb_dir, outputDir='templates'):
    """
    Maps PDB IDs to the main protein DF and copies the files locally.
    """
    # 1. Load the Excel mapping
    print('Reading template key file...')
    pdb_map_df = pd.read_excel(excel_path)
    
    # 2. Merge PDB IDs into the main DataFrame
    # We only need the UniProt_ID and the PDB structure file name (column titled "ID")
    pdb_map_df = pdb_map_df[['UniProt_ID', 'ID']].rename(columns={'UniProt_ID': 'UniProt ID','ID': 'PDB_File_Path'})
    print('Read template key file.')
    
    # Left join ensures we keep all proteins, even those without a PDB structure
    print('Mapping templates to proteins...')
    main_df = main_df.merge(pdb_map_df, left_on='UniProt ID', right_on='UniProt ID', how='left')
    
    found_count = main_df['PDB_File_Path'].notna().sum()
    print(f'Mapped {found_count} templates to proteins out of {len(main_df)} total proteins.')
    
    # 3. Handle File Copying
    if not os.path.exists(outputDir):
        print(f'Creating output directory for template files: "{outputDir}" ...')
        os.makedirs(outputDir)
        
    print(f"Copying PDB template files to {outputDir}...")
    
    def copy_if_exists(pdbFilePath: str):
        if pd.isna(pdbFilePath):
            return None
        
        # Files are named exactly after the 'ID' column + .pdb extension
        src_file = os.path.join(source_pdb_dir, f"{pdbFilePath}.pdb")
        uniProtID, pdbID = pdbFilePath.split('_')[1:3]
        pdbID = pdbID.upper()
        dst_file = os.path.join(outputDir, f"{uniProtID}_{pdbID}.pdb")
        
        if os.path.exists(src_file):
            shutil.copy(src_file, dst_file)
            return dst_file # Return the local path for the YAML generator
        else:
            return None

    # Apply the copy function and store the local path in a new column
    main_df['PDB_File_Path'] = main_df['PDB_File_Path'].apply(copy_if_exists)
    
    print(f"Successfully mapped and copied {found_count} PDB templates to '{outputDir}'.")
    
    return main_df

def main():
    # read protein file into df's
    proteins = pd.read_csv(PROTEIN_FILEPATH, index_col='ID')
    # trim unnecessary data from the protein df
    proteins = proteins[['Name', 'UniProt ID']].drop_duplicates()
    
    print(f'Read and processed {len(proteins)} proteins.')
    
    # get protein sequences by UniProt ID and add them as another column to the protein df
    proteins = fetch_uniprot_sequences(proteins)
    proteins = mapPDBtemplates(proteins, TEMPLATE_KEY_FILEPATH, TEMPLATES_DIR, os.path.join('output', 'templates'))
    proteins.to_csv(PROTEIN_OUTPUT_FILEPATH)
    
    print(f'Saved proteins with sequences and templates to {PROTEIN_OUTPUT_FILEPATH}')

if __name__ == '__main__':
    main()