import pandas as pd
import os
from itertools import product
import yaml

PROTEIN_SEQUENCE_FILEPATH = r'proteins_with_sequences.csv'
LIGAND_DIR = r'natural_products'
OUTPUT_DIR = r'output'

def filterDf(df: pd.DataFrame, columns: list[str]):
    '''
    Removes rows of a specified df containing empty strings or NaN values in the specified columns and returns the df with index reset.
    '''
    
    df = df.dropna(subset=columns)
    for column in columns:
        df = df[df[column].str.strip().astype(bool)]
        
    return df.reset_index(drop=True)

def generateBoltzYamls(proteinDf: pd.DataFrame, ligandDf: pd.DataFrame, outputDir, cap = float('inf')):
    """
    Generates Boltz YAML files for all combinations of protein and ligand given DataFrames containing the respective proteins and ligands.
    """
    
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    
    count = 0
    print(f"Generating Boltz YAML files...")

    # Cartesian product of cleaned DataFrames
    for (_, p_row), (_, l_row) in product(proteinDf.iterrows(), ligandDf.iterrows()):
        
        if count >= cap:
            break
            
        uniprot_id = p_row['UniProt ID']
        protein_seq = p_row['Sequence']
        np_mrd_id = l_row['NP_MRD_ID']
        smiles = l_row['SMILES']
        
        # Construct the exact structure you provided
        yaml_data = {
            "sequences": [
                {
                    "protein": {
                        "id": "A",
                        "sequence": protein_seq
                    }
                },
                {
                    "ligand": {
                        "id": "B",
                        "smiles": smiles
                    }
                }
            ],
            "properties": [
                {
                    "affinity": {
                        "binder": "B"
                    }
                }
            ]
        }
        
        # Filename based on indices to ensure uniqueness
        # e.g., P45059_L42.yaml
        filename = f"{uniprot_id}_{np_mrd_id}.yaml"
        file_path = os.path.join(outputDir, filename)
        
        # Write the YAML file
        with open(file_path, 'w') as f:
            # We use default_flow_style=False to keep it human-readable
            # and sort_keys=False to maintain the order (sequences -> properties)
            yaml.dump(yaml_data, f, default_flow_style=False, sort_keys=False)
            
        count += 1
        if count % 100 == 0:
            print(f"Progress: {count} files generated.")

    print(f"\nSuccess: {count} YAML files written to '{outputDir}/'")

def main():
    # read the proteins (w/ sequences) and ligands into df's
    print('Reading data...')
    proteins = pd.read_csv(PROTEIN_SEQUENCE_FILEPATH, index_col='ID')
    ligands = pd.concat(
        (pd.read_csv(f) for f in os.scandir(LIGAND_DIR)),
        ignore_index=True
        )
    
    print(f'Read {len(proteins)} proteins and {len(ligands)} ligands.')
    
    # filter the data
    # by default, removes any proteins with blank sequences or ligands with blank structures
    # can also implement functionality to only keep specific protein or ligand ids
    proteins = filterDf(proteins, ['UniProt ID', 'Sequence'])
    ligands = filterDf(ligands, ['SMILES'])
    
    print(f'After filtering: {len(proteins)} proteins and {len(ligands)} ligands remain.')
    
    # generate boltz-ready yaml files using the compiled protein and ligand data
    generateBoltzYamls(proteins, ligands, OUTPUT_DIR, cap=50)

if __name__ == '__main__':
    main()