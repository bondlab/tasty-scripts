#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Renames 'region_output.gff3' inside each accession-named folder to '{accession}_regionoutput.gff3',
then copies them into a destination directory (default: signal_<source_dirname>).
"""

import os
import shutil
import re
import argparse
from pathlib import Path

def rename_and_copy_files(source_dir, destination_dir):
    source_dir = Path(source_dir).resolve()

    # Create default destination if none provided
    if destination_dir is None:
        destination_dir = source_dir.parent / f"signal_{source_dir.name}"
    else:
        destination_dir = Path(destination_dir).resolve()

    destination_dir.mkdir(parents=True, exist_ok=True)
    print(f"📁 Copying renamed files to: {destination_dir}")

    for accession in os.listdir(source_dir):
        folder_path = source_dir / accession
        if not folder_path.is_dir():
            continue

        old_file = folder_path / "region_output.gff3"
        
        # Extract clean accession name from folder using regex
        match = re.search(r"(G[ACF]A?_\d+\.\d+)", accession)
        if match:
            accession_base = match.group(1)
        else:
            print(f"⚠️ Could not extract accession from folder name: {accession}")
            continue  # skip if no accession found
        
        new_file_name = f"{accession_base}_regionoutput.gff3"
        
        new_file_path = folder_path / new_file_name

        if old_file.is_file():
            os.rename(old_file, new_file_path)
            shutil.copy2(new_file_path, destination_dir / new_file_name)
            print(f"✅ Renamed and copied: {new_file_name}")
        else:
            print(f"⚠️ No region_output.gff3 found in {accession}")

def main():
    parser = argparse.ArgumentParser(description="Rename and copy region_output.gff3 files from subfolders.")
    parser.add_argument("source", help="Directory containing accession-named subfolders")
    parser.add_argument("-o", "--output", help="Destination directory for renamed files (default: signal_<source_foldername>)")
    args = parser.parse_args()

    rename_and_copy_files(args.source, args.output)

if __name__ == "__main__":
    main()