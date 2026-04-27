# bioin-401
Code for OT-Finder &amp; ND-Finder
# Welcome to OT-Finder

This repository contains the end-to-end OT-Finder workflow used to:
- prepare drug-target input data,
- generate Boltz YAMLs,
- run large batches through the cluster ("meta farm"),
- and store/export scored results.

The guide below combines the SEA pipeline and downstream result handling into one step-by-step runbook.

## 0) Start in the project folder

```bash
cd ~/ot-finder/bioin-401
```

## 1) Confirm your input CSV has required columns

Before running YAML generation, make sure your input file contains all of:
- `DrugBank ID`
- `Query Smiles`
- `UniProt ID`
- `Target Sequence`

## 2) If fields are split across files, build one unified CSV first

If your data is split (for example, one file with drug info and another with protein info), create a single combined CSV before YAML generation.

The existing script `yaml-generation-sea/build_sea_sample.py` is the reference implementation for this:
- it reads SEA-filtered rows,
- maps DrugBank IDs + SMILES from `yaml-generation/structure_links.csv`,
- and fetches UniProt accessions + sequences,
- then writes `yaml-generation-sea/sea_data.csv`.

Run it as-is for SEA-based inputs:

```bash
python yaml-generation-sea/build_sea_sample.py
```

Optional sample mode:

```bash
python yaml-generation-sea/build_sea_sample.py --sample
```

If you are using custom non-SEA split files, follow the same pattern as `build_sea_sample.py`: identify the drug source file, protein source file, merge on your chosen key(s), and write one CSV containing the four required columns above.

## 3) If importing fresh SEA output, clean/filter it first

When bringing in a raw SEA result file, run:

```bash
python yaml-generation-sea/process_sea_results.py
```

This performs important preprocessing before YAML generation:
- applies a `Max Tc >= 0.35` threshold,
- normalizes/filter mammalian homologs (via `filter_mammalian_targets.py`),
- keeps only human, bacterial, or viral targets,
- writes `yaml-generation-sea/sea_results_filtered.csv`.

After this, build/rebuild `sea_data.csv` (Step 2) so downstream YAML generation gets clean input.

## 4) Generate YAMLs

Run:

```bash
python yaml-generation-sea/yaml_generator_sea.py
```

This reads `yaml-generation-sea/sea_data.csv` and writes YAMLs to:
- `~/scratch/ot-finder/sea_yamls`

## 5) Run the meta farm (cluster batch jobs)
```bash
farm_init.run Boltz2_run
```

Submit your batch scripts (example shown for SEA):

```bash
sbatch run_boltz_sea.sh
```

If you split SEA into chunks, submit `run_boltz_sea1.sh`, `run_boltz_sea2.sh`, etc.

Related runners:
- `run_boltz_ot.sh` for known off-target YAML runs.
- `run_boltz_avoidome.sh` / `run_boltz_avoidome1.sh` for avoidome runs.

## 6) After meta farm completes, ingest + store results

Use `yaml-generation/manage_results2_storage.py` to store results in SQLite and export CSV:

```bash
python yaml-generation/manage_results2_storage.py \
  --yaml-dir yaml-generation/output/yamls \
  --results-dir yaml-generation/results2 \
  --db-path yaml-generation/boltz_results.sqlite \
  --out-csv yaml-generation/boltz_results_export.csv
```

This is the main script to put results into a database and generate the exported results CSV.

## 7) SEA-specific export helper (optional but useful)

For SEA result folders under `~/scratch/ot-finder/sea_results`, use:

```bash
python yaml-generation-sea/export_boltz_sea_results.py
```

This updates:
- `yaml-generation-sea/boltz_sea_results.sqlite`
- `yaml-generation-sea/boltz_sea_results.csv`

## 8) Additional scripts worth knowing

- `yaml-generation-sea/filter_mammalian_targets.py`: standalone mammalian ortholog normalization utility.

## Typical SEA end-to-end command order

```bash
cd ~/ot-finder/bioin-401
python yaml-generation-sea/process_sea_results.py
python yaml-generation-sea/build_sea_sample.py
python yaml-generation-sea/yaml_generator_sea.py
sbatch run_boltz_sea.sh
# wait for jobs to finish
python yaml-generation-sea/export_boltz_sea_results.py
```

