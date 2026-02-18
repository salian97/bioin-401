'''
(not complete)
Given a dictionary (TARGET_MAP) of UniProt IDs, associated holo and apo PDB structures, and CCD codes for the ligands, generates Boltz-2-ready .yaml files to test the effect of template conditioning on Boltz-2 prediction accuracy
'''

import os
import json
import urllib.request
import ssl

# SSL context to bypass certificate issues on clusters
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE


TARGET_MAP = {
    "P00519": ["2GQG", "1OPJ", "DIT"], "P00533": ["1M17", "1XKK", "IRE"],
    "P07550": ["3SN6", "2RH1", "EPN"], "P03372": ["1ERE", "3ERT", "EST"],
    "P56817": ["1FKN", "2OHM", "992"], "P37231": ["1FM6", "1ZGY", "BRL"],
    "P29274": ["2YDV", "3EML", "NEC"], "P21554": ["6NI2", "5TGZ", "FUB"],
    "P24941": ["1HCK", "3PXF", "ATP"], "Q16539": ["1A9U", "1KV2", "SB2"],
    "P06213": ["1IR3", "4IBM", "ANP"], "P00734": ["1PPB", "1H8D", "0IH"],
    "P00742": ["2W26", "1P0S", "RIV"], "P27487": ["1X70", "6I7U", "STG"],
    "P41143": ["6V3Z", "4EJ4", "DPI"], "P10275": ["1E3G", "1Z95", "DHT"],
    "P04150": ["1M2Z", "1NHZ", "DEX"], "P35968": ["4AG8", "4ASD", "AXI"],
    "P12931": ["1A07", "1QCF", "STU"], "P07949": ["2IVU", "4UXP", "V02"],
    "P04626": ["3PP0", "3RCD", "SYR"], "P11309": ["2O3P", "5N4O", "QUC"],
    "P19793": ["1FBY", "2PZN", "RET"], "P08913": ["6KUX", "1HOD", "RS7"],
    "P41597": ["7XRL", "5T1A", "CCL"], "P42574": ["1GFW", "3IBF", "DEV"],
    "P17948": ["1FLT", "3V2A", "FVR"], "Q15596": ["1M2Z", "1GWQ", "NCO"],
    "P62509": ["2P7G", "1S9Q", "GSK"], "P41145": ["6B73", "4DJH", "6Y0"]
}

def fetch_api(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    try:
        with urllib.request.urlopen(req, context=ctx, timeout=15) as r:
            return json.loads(r.read().decode())
    except Exception as e:
        return None

def get_config_v5(uniprot_id, pdb_id, lig_id):
    pdb_id = pdb_id.upper()
    # 1. Get Entry Info (contains list of entities)
    entry = fetch_api(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")
    if not entry: return None, None, None
    
    target_seq = None
    l_type, l_val = None, None
    
    # 2. Iterate through Polymer Entities to find Receptor and Peptide Ligands
    polymers = entry.get('rcsb_entry_container_identifiers', {}).get('polymer_entity_ids', [])
    for eid in polymers:
        p_data = fetch_api(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{eid}")
        if not p_data: continue
        
        up_ids = p_data.get('rcsb_polymer_entity_container_identifiers', {}).get('uniprot_ids', [])
        seq = p_data.get('entity_poly', {}).get('pdbx_seq_one_letter_code')
        
        # Check if this is our target
        if up_ids and uniprot_id in up_ids:
            target_seq = seq
        # Check if this is a peptide ligand (not our target, no other UniProt, and short)
        elif not up_ids and seq and len(seq) < 60 and not l_val:
            l_val = seq
            l_type = "peptide"

    # 3. Check for Small Molecule if no peptide ligand found
    if not l_type:
        chem = fetch_api(f"https://data.rcsb.org/rest/v1/core/chemcomp/{lig_id.upper()}")
        if chem:
            l_val = chem.get('rcsb_chem_comp_descriptor', {}).get('smiles')
            if l_val: l_type = "small_molecule"

    return target_seq, l_type, l_val

def write_yaml(path, seq, lt, lv, template_path=None, force=False):
    with open(path, 'w') as f:
        f.write("sequences:\n")
        f.write(f"  - protein:\n      id: A\n      sequence: \"{seq}\"\n")
        if lt == "small_molecule":
            f.write(f"  - ligand:\n      id: B\n      smiles: \"{lv}\"\n")
        elif lt == "peptide":
            f.write(f"  - protein:\n      id: B\n      sequence: \"{lv}\"\n")
        if template_path:
            f.write("templates:\n  - protein:\n      id: A\n")
            f.write(f"      path: \"../{template_path}\"\n")
            if force: f.write("      force: true\n      threshold: 1.0\n")

# Execution
for g in ["group_1", "group_2", "group_3"]: os.makedirs(g, exist_ok=True)
os.makedirs("templates", exist_ok=True)

print("Running REST-based Setup (v5)...")
for uniprot, ids in TARGET_MAP.items():
    print(f"Target: {uniprot} ({ids[0]})...", end=" ", flush=True)
    t_seq, lt, lv = get_config_v5(uniprot, ids[0], ids[2])
    
    if not t_seq or not lv:
        print("FAILED.")
        continue

    write_yaml(f"group_1/{uniprot}.yaml", t_seq, lt, lv)
    
    apo_pdb = ids[1].upper()
    apo_file = f"templates/{apo_pdb}.pdb"
    if not os.path.exists(apo_file):
        try: urllib.request.urlretrieve(f"https://files.rcsb.org/download/{apo_pdb}.pdb", apo_file)
        except: pass
        
    write_yaml(f"group_2/{uniprot}.yaml", t_seq, lt, lv, apo_file)
    write_yaml(f"group_3/{uniprot}.yaml", t_seq, lt, lv, apo_file, force=True)
    print("SUCCESS.")