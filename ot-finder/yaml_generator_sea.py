import argparse
import os

import pandas as pd
import yaml
from rdkit import Chem

# Local input tables.
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_FILEPATH = os.path.join(BASE_DIR, "sea_data.csv")
KNOWN_OFF_TARGETS_FILEPATH = os.path.join(BASE_DIR, "known_off_targets.csv")
AVOIDOME_FILEPATH = os.path.join(BASE_DIR, "avoidome_sample.csv")
# Output locations on scratch for each workflow mode.
OUTPUT_DIR = os.path.normpath(os.path.expanduser("~/scratch/ot-finder/sea_yamls"))
PILOT_RESULTS_DIR = os.path.normpath(os.path.expanduser("~/scratch/ot-finder/sea_results"))
OT_YAMLS_DIR = os.path.normpath(os.path.expanduser("~/scratch/ot-finder/ot_yamls"))
OT_RESULTS_DIR = os.path.normpath(os.path.expanduser("~/scratch/ot-finder/ot_results"))
AVOIDOME_YAMLS_DIR = os.path.normpath(os.path.expanduser("~/scratch/ot-finder/avoidome_yaml"))
AVOIDOME_RESULTS_DIR = os.path.normpath(os.path.expanduser("~/scratch/ot-finder/avoidome_results"))


def filter_df(df: pd.DataFrame, columns: list) -> pd.DataFrame:
    """Remove rows with missing or empty values in specified columns."""
    # Drop true NaNs first, then trim/guard against whitespace-only strings.
    df = df.dropna(subset=columns)
    for col in columns:
        df = df[df[col].astype(str).str.strip().astype(bool)]
    return df.reset_index(drop=True)


def save_csv(df: pd.DataFrame, path: str) -> None:
    df.to_csv(path, index=False)


def generate_boltz_yamls(
    protein_ligand_pairs: pd.DataFrame, output_dir: str, msa_results_dir: str
) -> None:
    """
    Write one YAML per row under output_dir. Dataframe must already be ordered by
    UniProt ID so pilot MSA reuse matches consecutive rows per protein.
    First YAML for each UniProt is pilot (no msa); later ones point to that pilot's MSA.
    """
    os.makedirs(output_dir, exist_ok=True)
    print("Generating Boltz YAML files...")
    # Tracks the first YAML stem seen for each UniProt (the "pilot" run).
    pilot_stem_by_uniprot: dict[str, str] = {}
    # Also used as the file index prefix in generated YAML names.
    count = 0

    for _, row in protein_ligand_pairs.iterrows():
        # Normalize to strings so CSV mixed types do not break YAML emission.
        uniprot_id = str(row["UniProt ID"]).strip()
        drugbank_id = str(row["DrugBank ID"]).strip()
        smiles = str(row["Query Smiles"]).strip()
        protein_seq = str(row["Target Sequence"]).strip()

        # Validate the ligand; invalid SMILES are skipped rather than failing the run.
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Skipping {drugbank_id}/{uniprot_id}: invalid SMILES")
            continue

        atom_count = mol.GetNumAtoms()

        yaml_data = {
            "sequences": [
                {"protein": {"id": "A", "sequence": protein_seq}},
                {"ligand": {"id": "B", "smiles": smiles}},
            ]
        }

        # Boltz affinity property is only emitted for ligands up to 128 atoms.
        if atom_count <= 128:
            yaml_data["properties"] = [{"affinity": {"binder": "B"}}]
        else:
            print(f"{drugbank_id}: {atom_count} atoms → affinity skipped")

        # File naming convention is consumed by downstream ingestion/export scripts.
        filename = f"{count}_{uniprot_id}_{drugbank_id}.yaml"
        file_path = os.path.join(output_dir, filename)
        file_stem = os.path.splitext(filename)[0]

        if uniprot_id in pilot_stem_by_uniprot:
            pilot_stem = pilot_stem_by_uniprot[uniprot_id]
            # Non-pilot rows for this protein reuse the pilot MSA CSV path.
            msa_path = os.path.join(
                msa_results_dir,
                f"boltz_results_{pilot_stem}",
                "msa",
                f"{pilot_stem}_0.csv",
            )
            yaml_data["sequences"][0]["protein"]["msa"] = msa_path
        else:
            # First time this UniProt appears: keep it as pilot (no msa field).
            pilot_stem_by_uniprot[uniprot_id] = file_stem

        with open(file_path, "w") as f:
            yaml.dump(yaml_data, f, default_flow_style=False, sort_keys=False)

        count += 1
        if count % 100 == 0:
            print(f"Progress: {count} files generated.")

    print(f"\nSuccess: {count} YAML files written to {output_dir}")
    print(f"MSA results base directory = {msa_results_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate Boltz YAMLs from SEA or known off-target CSVs.")
    parser.add_argument(
        "--known-off-targets",
        action="store_true",
        help=f"Read {os.path.basename(KNOWN_OFF_TARGETS_FILEPATH)}, write YAMLs to {OT_YAMLS_DIR}, "
        f"MSA paths under {OT_RESULTS_DIR}. Does not modify sea_data.csv.",
    )
    parser.add_argument(
        "--avoidome",
        action="store_true",
        help=f"Read {os.path.basename(AVOIDOME_FILEPATH)}, write YAMLs to {AVOIDOME_YAMLS_DIR}, "
        f"MSA paths under {AVOIDOME_RESULTS_DIR}. Does not modify the CSV.",
    )
    args = parser.parse_args()

    if args.known_off_targets and args.avoidome:
        raise SystemExit("Use only one of --known-off-targets or --avoidome.")

    # Mode switch controls both input schema and output destinations.
    if args.known_off_targets:
        input_path = KNOWN_OFF_TARGETS_FILEPATH
        output_dir = OT_YAMLS_DIR
        msa_results_dir = OT_RESULTS_DIR
        rewrite_input_csv = False
    elif args.avoidome:
        input_path = AVOIDOME_FILEPATH
        output_dir = AVOIDOME_YAMLS_DIR
        msa_results_dir = AVOIDOME_RESULTS_DIR
        rewrite_input_csv = False
    else:
        input_path = INPUT_FILEPATH
        output_dir = OUTPUT_DIR
        msa_results_dir = PILOT_RESULTS_DIR
        rewrite_input_csv = True

    print(f"Reading {input_path}...")
    df = pd.read_csv(input_path)

    if args.known_off_targets or args.avoidome:
        # Normalize alternate input headers to the SEA names expected downstream.
        rename = {
            "Drug SMILES": "Query Smiles",
            "Protein sequence": "Target Sequence",
        }
        df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})

    required = ["DrugBank ID", "UniProt ID", "Query Smiles", "Target Sequence"]
    missing_cols = [c for c in required if c not in df.columns]
    if missing_cols:
        raise ValueError(f"Input CSV missing columns: {missing_cols}")

    # Stable sort by protein first so pilot MSA reuse is deterministic.
    sort_by_uniprot = ["UniProt ID", "DrugBank ID"]
    if "Target ID" in df.columns:
        sort_by_uniprot.append("Target ID")
    if "Drug Name" in df.columns:
        sort_by_uniprot.append("Drug Name")

    df = df.sort_values(by=sort_by_uniprot, na_position="last", kind="mergesort")
    df = df.reset_index(drop=True)

    if rewrite_input_csv:
        # SEA mode rewrites the input table in UniProt order for reproducibility.
        print(f"Writing {input_path} ordered by UniProt ID ({len(df)} rows)...")
        save_csv(df, input_path)

    # Final safety filter before YAML emission.
    working = filter_df(df, required)
    print(f"{len(working)} protein-ligand pairs after filtering")

    generate_boltz_yamls(working, output_dir, msa_results_dir)

    if rewrite_input_csv:
        # Restore DrugBank-first ordering for easier review/lookup after generation.
        sort_by_drugbank = ["DrugBank ID", "UniProt ID"]
        if "Target ID" in df.columns:
            sort_by_drugbank.append("Target ID")
        df_out = df.sort_values(by=sort_by_drugbank, na_position="last", kind="mergesort").reset_index(
            drop=True
        )
        print(f"Writing {input_path} ordered by DrugBank ID ({len(df_out)} rows)...")
        save_csv(df_out, input_path)


if __name__ == "__main__":
    main()
