import pandas as pd
import os
import glob
import re
from pathlib import Path

SEA_RESULTS_DIR = r'sea_results'
SPECLIST_PATH = r'speclist.txt'
OUTPUT_DIR = r'outputs'

def parse_uniprot_speclist(file_path):
    data = []
    current_entry = None
    
    # Regex patterns for the three types of names
    # (.*?) is a non-greedy match to stop before the next tag
    scientificNamePattern = re.compile(r"N=(.*?)(?:[CS]=|$)")
    commonNamePattern = re.compile(r"C=(.*?)(?:[NS]=|$)")
    synonymPattern = re.compile(r"S=(.*?)(?:[NC]=|$)")

    with open(file_path, 'r', encoding='utf-8') as f:
        start_parsing = False
        for line in f:
            if '_____' in line:
                start_parsing = True
                continue
            
            if not start_parsing or not line.strip():
                continue
            
            if line.startswith('====='):
                break

            # If the first column (0-5) has text, it's a NEW species entry
            if line[0:5].strip() != "":
                # Save the finished previous entry before starting a new one
                if current_entry:
                    data.append(process_entry(current_entry, scientificNamePattern, commonNamePattern, synonymPattern))
                
                # Initialize new entry dictionary
                current_entry = {
                    'SpeciesCode': line[0:5].strip(),
                    'Kingdom': line[6:7].strip(),
                    'TaxonNode': int(line[8:15].strip()),
                    'RawNames': line[16:].strip()
                }
            else:
                # If it starts with spaces, append to the current entry's names
                if current_entry:
                    current_entry['RawNames'] += " " + line.strip()

        # Append the very last entry in the file
        if current_entry:
            data.append(process_entry(current_entry, scientificNamePattern, commonNamePattern, synonymPattern))

    return pd.DataFrame(data).set_index('SpeciesCode')

def process_entry(entry: dict, sci_p: re.Pattern[str], com_p: re.Pattern[str], syn_p: re.Pattern[str]):
    """Helper to clean strings and extract N, C, and S components."""
    
    rawNames = entry['RawNames']
    
    # Extract using regex
    sci_match = sci_p.search(rawNames)
    com_match = com_p.search(rawNames)
    syn_match = syn_p.search(rawNames)
    
    return {
        'SpeciesCode': entry['SpeciesCode'],
        'Kingdom': entry['Kingdom'],
        'TaxonNode': int(entry['TaxonNode']),
        'ScientificName': sci_match.group(1).strip() if sci_match else None,
        'CommonName': com_match.group(1).strip() if com_match else None,
        'Synonyms': syn_match.group(1).strip() if syn_match else None
    }
    
def filterSpecies(seaResults: pd.DataFrame, speciesList: pd.DataFrame):
    # 1. Extract SpeciesCode from Target ID (e.g., '5HT2C_HUMAN' -> 'HUMAN')
    # We use .str.split('_').str[-1] to get the part after the last underscore
    seaResults['SpeciesCode'] = seaResults['Target ID'].str.split('_').str[-1]
    
    # 2. Map Kingdom from the species list
    # Since SpeciesCode is the index of df_species, .map() is very efficient
    seaResults['Kingdom'] = seaResults['SpeciesCode'].map(speciesList['Kingdom'])
    
    # 3. Filter for Human, Bacterial, or Viral targets
    # Note: 'HUMAN' is in Kingdom 'E' (Eukaryota), so we check for it specifically
    # to avoid including other Eukaryotes like RAT or MOUSE.
    speciesFilter = (
        (seaResults['SpeciesCode'] == 'HUMAN') | 
        (seaResults['Kingdom'] == 'B') | 
        (seaResults['Kingdom'] == 'V')
    )
    
    filteredSeaResults = seaResults[speciesFilter].copy()
    
    # Optional: Clean up the helper columns if you don't need them in the final output
    return filteredSeaResults.drop(columns=['SpeciesCode', 'Kingdom'])

def main():
    print('Reading SEA search results...')
    # build df with all the SEA result files
    allCSVs = glob.glob(os.path.join(SEA_RESULTS_DIR, "*"))
    seaResults = pd.concat([pd.read_csv(f) for f in allCSVs], ignore_index=True)
    print(f'Read {len(seaResults)} rows of SEA results.')
    
    # save a list of all of the original natural products
    ''' get the list that Jordan will put together and compare it to the SMILES codes in the SEA results '''
    
    print('Filtering out proteins with maxTC < 0.5 ...')
    seaResults = seaResults[seaResults['Max Tc'] >= 0.5]
    print(f'After filtering, {len(seaResults)} rows remain.')
    
    print('Loading UniProt species vocabulary data ...')
    speciesList = parse_uniprot_speclist(SPECLIST_PATH)
    print(f'UniProt species vocabulary data loaded. ({len(speciesList)} rows)')
    
    print(seaResults)
    print(speciesList)
    
    print('Filtering for only human, bacterial, or viral targets ...')
    seaResults = filterSpecies(seaResults, speciesList).reset_index()
    print(f'Filtered for only human, bacterial, or viral targets. ({len(seaResults)} rows)')
    
    print('Saving filtered SEA results ...')
    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
    outputPath = Path(OUTPUT_DIR).joinpath('processed_sea_results.csv')
    seaResults.to_csv(outputPath)
    print(f'Filtered SEA results saved to "{outputPath}" !')

if __name__ == '__main__':
    main()