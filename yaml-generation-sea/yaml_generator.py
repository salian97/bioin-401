import pandas as pd
import os
import yaml

OUTPUT_DIR = r'outputs'
INPUT_FILEPATH = os.path.join(OUTPUT_DIR, 'processed_sea_results_with_uniprot_data.csv')
OUTPUT_FILEPATH = os.path.join(OUTPUT_DIR, 'yamls')

# def filterDf(df: pd.DataFrame, columns: list[str]):
#     '''
#     Removes rows of a specified df containing empty strings or NaN values in the specified columns and returns the df with index reset.
#     '''
    
#     df = df.dropna(subset=columns)
#     for column in columns:
#         df = df[df[column].str.strip().astype(bool)]
        
#     return df.reset_index(drop=True)

def generateBoltzYamls(proteinLigandPairs: pd.DataFrame, outputDir, cap=float('inf'), chunk_size=float('inf')):
    """
    Generates Boltz YAML files with sub-directory chunking and proper relative paths.
    """
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    
    count = 0
    print(f"Generating Boltz YAML files...")

    # Determine if we should use subdirectories
    use_chunks = cap > chunk_size

    for _, row in proteinLigandPairs.iterrows():
        if count >= cap:
            break
            
        uniprot_id = row['UniProt ID']
        protein_seq = row['Sequence']
        template_path = row['PDB_File_Path']
        compoundId = row['Query ID']
        smiles = row['Query Smiles']
        
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
            
            filename = f"{count}_{uniprot_id}_{compoundId}_T.yaml"
        else:
            filename = f"{count}_{uniprot_id}_{compoundId}.yaml"
            
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
    # read processed protein-ligand pairs into a df
    print('Reading data...')
    proteinLigandPairs = pd.read_csv(INPUT_FILEPATH)
    
    print(f'Read {len(proteinLigandPairs)} protein-ligand pairs.')
    
    # generate boltz-ready yaml files using the compiled protein and ligand data
    generateBoltzYamls(proteinLigandPairs, os.path.join(OUTPUT_DIR, 'yamls'), cap=100, chunk_size=5)

if __name__ == '__main__':
    main()