'''
(not complete)
Given a UniProt ID, identifies the "best" experimental PDB structure available and generates a YAML file.
In the program's current state, it identifies and parses the different PDB structures available for a given UniProt ID. To complete the program, the next step would be to implement a scoring method to select the "best" PDB structure, and then fetch that PDB structure using the RCSB PDB API.
'''

import sys
import json
import requests
import pandas as pd

UNIPROT_API = 'https://rest.uniprot.org/uniprotkb/'

# helper function - wraps around the requests library and adds error handling
def getUrl(url, **kwargs):
    response = requests.get(url, **kwargs) # send a GET request to the given URL
    
    # error handling
    if not response.ok:
        print(response.text) # print error message to the terminal
        response.raise_for_status() # raises an http error
        sys.exit() # exits the program
    
    return response # if the response was successful, return it

def pdbXrefsToDf(pdbXrefs):
    # initialize an empty list to add all the data to, which will then be converted to a df
    rows = []
    
    # iterate through each pdb cross-reference in the json
    for xref in pdbXrefs:
        # initialize the entry by adding its pdb id to a new dict
        entry = {'PDB_ID': xref['id']}
        
        
        for prop in xref.get('properties', []):
            entry[prop["key"]] = prop["value"]
            
        rows.append(entry)
    
    df = pd.DataFrame(rows)
    
    # remove the resolution units, NaN if no resolution for that xref
    df['Resolution'] = pd.to_numeric(df['Resolution'].str.replace(' A', ''), errors='coerce')
    
    # split the chains and start/end residues into separate columns
    df[['Chains', 'Start', 'End']] = df['Chains'].str.extract(r'(.+)=(\d+)-(\d+)')
    df[['Start', 'End']] = df[['Start', 'End']].astype(int) # convert start and end residues to ints
    df['Chains'] = df['Chains'].str.split('/') # split the chains into a list
    
    return df

def domainsToDf(features):
    rows = []
    
    for feature in features:
        entry = {
            'Description': feature.get('description', ''),
            'Start': feature.get('location', {}).get('start', {}).get('value', ''),
            'End': feature.get('location', {}).get('end', {}).get('value', '')
        }
        rows.append(entry)
    
    return pd.DataFrame(rows)

def rankPDBs(pdbXrefs: pd.DataFrame, domains: pd.DataFrame):
    # 1. filter for only X-ray & EM structures
    filterCondition = pdbXrefs['Method'].isin(['X-ray', 'EM'])
    pdbXrefs = pdbXrefs[filterCondition].copy()
    
    # 2. calculate coverage of functionally-relevant protein domains for each pdb in a new column
    pdbXrefs['Coverage'] = pdbXrefs.apply(
        lambda row: calculate_pdb_coverage(row['Start'], row['End'], domains),
        axis=1
    )
    
    return pdbXrefs
    
def calculate_pdb_coverage(pdb_start, pdb_end, domains_df: pd.DataFrame):
    total_covered_residues = 0
    
    for _, domain in domains_df.iterrows():
        # Find the overlap between the PDB range and the current Domain range
        # Formula: Overlap Start is max of starts, Overlap End is min of ends
        overlap_start = max(pdb_start, domain['Start'])
        overlap_end = min(pdb_end, domain['End'])
        
        # If they overlap, the difference is positive
        if overlap_start <= overlap_end:
            total_covered_residues += (overlap_end - overlap_start + 1)
            
    return int(total_covered_residues)

def main():
    # xref_pdb gets all pdb xrefs, ft_domain gets a feature table of the protein domains
    r = getUrl(f"{UNIPROT_API}P00533?fields=xref_pdb,ft_domain")
    pdbXrefs = pdbXrefsToDf(r.json()['uniProtKBCrossReferences'])
    domains = domainsToDf(r.json()['features'])
        
    print(domains)
    
    print(pdbXrefs)
    
    pdbXrefs = rankPDBs(pdbXrefs, domains)
    
    print(pdbXrefs)
    
    # print(pdbXrefs)
    
    # filter for only x-ray structures
    
    # print(pdbXrefs)
    
    # testDf = pd.json_normalize(responseJson)
    # print(testDf)
    # print(json.dumps(r.json(), indent=2)[:1000])

if __name__ == '__main__':
    main()