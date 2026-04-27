import pandas as pd
import re
from pathlib import Path

from filter_mammalian_targets import filter_sea_results as apply_mammalian_target_filter

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_CSV = SCRIPT_DIR / "sea_results.csv"
SPECLIST_PATH = SCRIPT_DIR / "speclist.txt"
OUTPUT_CSV = SCRIPT_DIR / "sea_results_filtered.csv"
# Keep rows with Max Tc at or above this cutoff (drop lower).
MAX_TC_MIN = 0.35

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
    print(f"Reading SEA search results from {INPUT_CSV} ...")
    if not INPUT_CSV.is_file():
        raise FileNotFoundError(f"Input CSV not found: {INPUT_CSV}")
    seaResults = pd.read_csv(INPUT_CSV, low_memory=False)
    print(f"Read {len(seaResults)} rows.")

    print(f"Filtering out proteins with maxTC < {MAX_TC_MIN} ...")
    seaResults = seaResults[seaResults["Max Tc"] >= MAX_TC_MIN].copy()
    print(f"After Max Tc filter, {len(seaResults)} rows remain.")

    # Before kingdom filtering: collapse mammal orthologs to _HUMAN when no human row exists
    # for the same (Query ID, gene prefix), per filter_mammalian_targets.py.
    print("Normalizing mammalian ortholog targets (RAT/MOUSE/BOVIN/...) ...")
    fieldnames = list(seaResults.columns)
    row_dicts = seaResults.to_dict("records")
    row_dicts = apply_mammalian_target_filter(row_dicts, fieldnames)
    seaResults = pd.DataFrame(row_dicts, columns=fieldnames)
    print(f"After mammalian normalization, {len(seaResults)} rows.")

    print("Loading UniProt species vocabulary data ...")
    speciesList = parse_uniprot_speclist(str(SPECLIST_PATH))
    print(f"UniProt species vocabulary data loaded. ({len(speciesList)} rows)")

    print("Filtering for only human, bacterial, or viral targets ...")
    seaResults = filterSpecies(seaResults, speciesList).reset_index(drop=True)
    print(f"Filtered for only human, bacterial, or viral targets. ({len(seaResults)} rows)")

    print(f"Saving filtered SEA results to {OUTPUT_CSV} ...")
    seaResults.to_csv(OUTPUT_CSV, index=False)
    print(f'Filtered SEA results saved to "{OUTPUT_CSV}".')

if __name__ == '__main__':
    main()
