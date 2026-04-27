import yaml
import os

# --- CONFIGURATION: Define your source and target pathways ---

# The original Nibi-Wong components to be replaced
OLD_PREFIX = "/project/6002707/Serag"
OLD_FOLDER = "bioin401farm"

# The target components
NEW_PREFIX = "/project/def-gksw/elgamal"
NEW_FOLDER = "bioin401-fir-gksw"

# Using os.path.join for cross-platform compatibility (Windows/Linux)
SOURCE_DIR = os.path.join('filepath_transformer', '0_before')
TARGET_DIR = os.path.join('filepath_transformer', '1_after')

# -------------------------------------------------------------

def update_path(path):
    """Applies the string replacement to a specific path string."""
    if isinstance(path, str):
        new_path = path.replace(OLD_PREFIX, NEW_PREFIX)
        new_path = new_path.replace(OLD_FOLDER, NEW_FOLDER)
        return new_path
    return path

def recursive_path_fix(obj):
    """Recursively crawls the YAML dictionary to find 'msa', 'pdb', and 'cif' keys."""
    if isinstance(obj, dict):
        for key, value in obj.items():
            if key in ['msa', 'pdb', 'cif'] and isinstance(value, str):
                obj[key] = update_path(value)
            else:
                recursive_path_fix(value)
    elif isinstance(obj, list):
        for item in obj:
            recursive_path_fix(item)

def process_yamls():
    # We'll keep a counter just for your terminal output
    file_count = 0

    print(f"Starting transformation from {SOURCE_DIR} to {TARGET_DIR}...")

    # os.walk "crawls" through every subdirectory
    for root, dirs, files in os.walk(SOURCE_DIR):
        for filename in files:
            if filename.endswith('.yaml') or filename.endswith('.yml'):
                file_count += 1
                
                # 1. Get the full path to the source file
                source_path = os.path.join(root, filename)
                
                # 2. Determine the relative path (e.g., "001438/file.yaml")
                # This ensures we keep the directory tree the exact same
                rel_path = os.path.relpath(source_path, SOURCE_DIR)
                target_path = os.path.join(TARGET_DIR, rel_path)
                
                # 3. Create the target subdirectory if it doesn't exist
                target_subdir = os.path.dirname(target_path)
                if not os.path.exists(target_subdir):
                    os.makedirs(target_subdir)

                # 4. Process the file
                with open(source_path, 'r') as f:
                    try:
                        data = yaml.safe_load(f)
                        recursive_path_fix(data)

                        with open(target_path, 'w') as out_f:
                            yaml.safe_dump(data, out_f, default_flow_style=False, sort_keys=False)
                    
                    except yaml.YAMLError as exc:
                        print(f"Error processing {source_path}: {exc}")

    if file_count == 0:
        print(f"Still found 0 files. Double check that your YAMLs are inside: {os.path.abspath(SOURCE_DIR)}")
    else:
        print(f"Success! Transformed {file_count} files while maintaining the directory structure.")

if __name__ == "__main__":
    process_yamls()