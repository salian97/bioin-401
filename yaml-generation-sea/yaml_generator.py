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
    
    subDirKey = []
    globalFileCount = 0
    folderCounter = 0
    
    print(f"Generating Boltz YAML files...")

    # group the df by protein ID
    grouped = proteinLigandPairs.groupby('UniProt ID')
    
    for uniprot_id, group in grouped:
        if globalFileCount >= cap:
            break
        
        # slice the group into chunks based on chunk_size
        num_ligands = len(group)
        for i in range(0, num_ligands, chunk_size):
            if globalFileCount >= cap:
                break
            
            # create the zero-padded subdirectory name (0000, 0001, 0002, ...)
            folderName = str(folderCounter).zfill(len(str(len(proteinLigandPairs))))
            current_save_dir = os.path.join(outputDir, folderName)
            
            if not os.path.exists(current_save_dir):
                os.makedirs(current_save_dir)
            
            # add subdirectory to subdirectory key
            subDirKey.append({
                'subdir_code': folderName,
                'uniprot_id': uniprot_id
            })
            
            # get slice of ligands for this specific chunk
            chunk_slice = group.iloc[i : i+chunk_size]
            
            for _, row in chunk_slice.iterrows():
                if globalFileCount >= cap:
                    break
                
                protein_seq = row['Sequence']
                template_path = row['PDB_File_Path']
                compoundId = row['Query ID']
                smiles = row['Query Smiles']
                
                # 3. Build Boltz-2 YAML Structure
                yaml_data = {
                    "sequences": [
                        {"protein": {"id": "A", "sequence": protein_seq}},
                        {"ligand": {"id": "B", "smiles": smiles}}
                    ]
                }
                
                # Handling templates with your specific absolute path logic
                has_template = pd.notna(template_path) and os.path.exists(str(template_path))
                if has_template:
                    # Maintain the absolute path format for the cluster
                    yaml_template_path = '/project/6002707/Serag/envs/boltz/bioin401/bioin401farm/prod1/templates/' + os.path.basename(template_path)
                    ext = os.path.splitext(template_path)[1].lower()
                    template_key = ext[1:] # e.g., 'pdb' or 'cif'
                    
                    yaml_data["templates"] = [{template_key: yaml_template_path}]
                    filename = f"{globalFileCount}_{uniprot_id}_{compoundId}_T.yaml"
                else:
                    filename = f"{globalFileCount}_{uniprot_id}_{compoundId}.yaml"
                    
                yaml_data["properties"] = [{"affinity": {"binder": "B"}}]
                
                # Final write
                file_path = os.path.join(current_save_dir, filename)
                with open(file_path, 'w') as f:
                    yaml.dump(yaml_data, f, default_flow_style=False, sort_keys=False)
                
                globalFileCount += 1
                if globalFileCount % 100 == 0:
                    print(f"Progress: {globalFileCount} files generated.")
            
            # Increment folder counter after finishing a chunk
            folderCounter += 1
            
    # 4. Generate the Subdirectory Key CSV
    key_df = pd.DataFrame(subDirKey)
    key_path = os.path.join(OUTPUT_DIR, 'subdirectory_key.csv')
    key_df.to_csv(key_path, index=False)
    
    print(f"\nSuccess: {globalFileCount} YAMLs written across {folderCounter} folders.")
    print(f"Subdirectory key saved to: {key_path}")

def main():
    # read processed protein-ligand pairs into a df
    print('Reading data...')
    proteinLigandPairs = pd.read_csv(INPUT_FILEPATH)
    
    print(f'Read {len(proteinLigandPairs)} protein-ligand pairs.')
    
    # generate boltz-ready yaml files using the compiled protein and ligand data
    generateBoltzYamls(proteinLigandPairs, os.path.join(OUTPUT_DIR, 'yamls'), chunk_size=10)

if __name__ == '__main__':
    main()