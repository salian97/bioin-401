import pandas as pd
import requests
import io
import time

PROTEIN_FILEPATH = r'drugbank-targets.csv'
UNIPROT_API = r'https://rest.uniprot.org/uniprotkb/accessions'
OUTPUT_FILEPATH = r'proteins_with_sequences.csv'

def fetch_uniprot_sequences(df: pd.DataFrame, id_column='UniProt ID', chunk_size=500):
    """
    Fetches protein sequences from UniProt in batches and maps them to the DataFrame.
    """
    
    # 1. Get unique IDs to avoid redundant API calls
    unique_ids = df[id_column].unique().tolist()
    sequence_map = {}
    
    print(f"Starting retrieval for {len(unique_ids)} unique IDs...")

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

def main():
    # read protein and ligand csv files into df's
    proteins = pd.read_csv(PROTEIN_FILEPATH, index_col='ID')
    # trim unnecessary data from the protein df
    proteins = proteins[['Name', 'UniProt ID']]
    
    print(f'Read and processed {len(proteins)} proteins.')
    
    # get protein sequences by UniProt ID and add them as another column to the protein df
    proteins = fetch_uniprot_sequences(proteins)
    proteins.to_csv(OUTPUT_FILEPATH)
    
    print(f'Saved proteins with sequences to {OUTPUT_FILEPATH}')

if __name__ == '__main__':
    main()