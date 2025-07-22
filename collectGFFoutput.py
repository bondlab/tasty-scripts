#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Automatically detects file type in subfolders (region_output.gff3 or TMRs.gff3),
then renames and copies to a new directory:
- signal_<input_dir> for region_output.gff3
- beta_<input_dir> for TMRs.gff3
"""

import os
import shutil
import re
import argparse
from pathlib import Path

def detect_mode(source_dir):
    """Look in the first subdirectory to determine file mode."""
    for item in sorted(os.listdir(source_dir)):
        folder_path = source_dir / item
        if folder_path.is_dir():
            region_file = folder_path / "region_output.gff3"
            tmrs_file = folder_path / "TMRs.gff3"
            if region_file.is_file():
                return "signal"
            elif tmrs_file.is_file():
                return "beta"
            else:
                continue
    return None

def rename_and_copy_files(source_dir, destination_dir, mode):
    source_dir = Path(source_dir).resolve()

    # Decide filename and output dir prefix based on mode
    filename = "region_output.gff3" if mode == "signal" else "TMRs.gff3"
    if destination_dir is None:
        destination_dir = source_dir.parent / f"{mode}_{source_dir.name}"
    else:
        destination_dir = Path(destination_dir).resolve()

    destination_dir.mkdir(parents=True, exist_ok=True)
    print(f"📁 Copying renamed '{filename}' files to: {destination_dir}")

    for accession in os.listdir(source_dir):
        folder_path = source_dir / accession
        if not folder_path.is_dir():
            continue

        old_file = folder_path / filename

        # Extract accession (e.g., GCA_028697885.1)
        match = re.match(r"(G[ACF]F?_\d+\.\d+)", accession)
        if not match:
            match = re.search(r"(G[ACF]F?_\d+\.\d+)", accession)
        if match:
            accession_base = match.group(1)
        else:
            print(f"⚠️ Could not extract accession from folder name: {accession}")
            continue

        new_file_name = f"{accession_base}_{filename}"
        new_file_path = folder_path / new_file_name

        if old_file.is_file():
            os.rename(old_file, new_file_path)
            shutil.copy2(new_file_path, destination_dir / new_file_name)
            print(f"✅ Renamed and copied: {new_file_name}")
        else:
            print(f"⚠️ No {filename} found in {accession}")

def main():
    parser = argparse.ArgumentParser(description="Rename and copy GFF3 files from subfolders based on detected mode.")
    parser.add_argument("source", help="Directory containing accession-named subfolders")
    parser.add_argument("-o", "--output", help="Destination directory (default: signal_... or beta_...)")
    args = parser.parse_args()

    source_dir = Path(args.source).resolve()
    mode = detect_mode(source_dir)

    if not mode:
        print("❌ Could not detect mode (no region_output.gff3 or TMRs.gff3 found in any subfolder).")
        return

    print(f"🔍 Detected mode: {mode}")
    rename_and_copy_files(source_dir, args.output, mode)

if __name__ == "__main__":
    main()