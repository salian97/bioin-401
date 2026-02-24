import pandas as pd
import os
from itertools import product
import yaml

PROTEIN_SEQUENCE_FILEPATH = r'output_enzyme_gpcr_targets_with_sequences.csv'
LIGAND_FILEPATH = r'cluster_representatives_v1.csv'
OUTPUT_DIR = r'output'

def filterDf(df: pd.DataFrame, columns: list[str]):
    '''
    Removes rows of a specified df containing empty strings or NaN values in the specified columns and returns the df with index reset.
    '''
    
    df = df.dropna(subset=columns)
    for column in columns:
        df = df[df[column].str.strip().astype(bool)]
        
    return df.reset_index(drop=True)

def generateBoltzYamls(proteinDf: pd.DataFrame, ligandDf: pd.DataFrame, outputDir, cap=float('inf'), chunk_size=float('inf')):
    """
    Generates Boltz YAML files with sub-directory chunking and proper relative paths.
    """
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    
    count = 0
    print(f"Generating Boltz YAML files...")

    # Determine if we should use subdirectories
    use_chunks = cap > chunk_size

    for (_, p_row), (_, l_row) in product(proteinDf.iterrows(), ligandDf.iterrows()):
        if count >= cap:
            break
            
        uniprot_id = p_row['UniProt ID']
        protein_seq = p_row['Sequence']
        template_path = p_row['PDB_File_Path']
        np_mrd_id = l_row['NP_MRD_ID']
        smiles = l_row['SMILES']
        
        # 1. Determine current save directory and filename
        current_save_dir = outputDir
        if use_chunks:
            chunk_idx = str(count // chunk_size)
            current_save_dir = os.path.join(outputDir, chunk_idx)
            if not os.path.exists(current_save_dir):
                os.makedirs(current_save_dir)

        # 2. Start with the sequences block
        yaml_data = {
            "sequences": [
                {"protein": {"id": "A", "sequence": protein_seq}},
                {"ligand": {"id": "B", "smiles": smiles}}
            ]
        }
        
        # 3. Conditionally add templates with correct relative paths
        has_template = pd.notna(template_path) and os.path.exists(str(template_path))
        if has_template:
            # We calculate relative path from the CURRENT directory of the YAML
            rel_template_path = os.path.relpath(template_path, start=current_save_dir)
            yaml_data["templates"] = [{"pdb": rel_template_path}]
            
            filename = f"{count}_{uniprot_id}_{np_mrd_id}_T.yaml"
        else:
            filename = f"{count}_{uniprot_id}_{np_mrd_id}.yaml"
            
        # 4. Add properties
        yaml_data["properties"] = [{"affinity": {"binder": "B"}}]
        
        # Final write
        file_path = os.path.join(current_save_dir, filename)
        with open(file_path, 'w') as f:
            yaml.dump(yaml_data, f, default_flow_style=False, sort_keys=False)
            
        count += 1
        if count % 100 == 0:
            print(f"Progress: {count} files generated.")

    print(f"\nSuccess: {count} YAML files written to '{outputDir}'")

def main():
    # read the proteins (w/ sequences) and ligands into df's
    print('Reading data...')
    proteins = pd.read_csv(PROTEIN_SEQUENCE_FILEPATH)
    ligands = pd.read_csv(LIGAND_FILEPATH).drop_duplicates()
    
    print(f'Read {len(proteins)} proteins and {len(ligands)} ligands.')
    
    # filter the data
    # by default, removes any proteins with blank sequences or ligands with blank structures
    # can also implement functionality to only keep specific protein or ligand ids
    proteins = filterDf(proteins, ['UniProt ID', 'Sequence'])
    ligands = filterDf(ligands, ['SMILES'])
    
    print(f'After filtering: {len(proteins)} proteins and {len(ligands)} ligands remain.')
    
    # generate boltz-ready yaml files using the compiled protein and ligand data
    generateBoltzYamls(proteins, ligands, os.path.join(OUTPUT_DIR, 'yamls'), cap=10000, chunk_size=1000)

if __name__ == '__main__':
    main()