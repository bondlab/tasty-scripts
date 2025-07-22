#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 6 16:40:13 2025

@author: daniel

Give it a directory.

It looks for folders, and inside each folder, it will find all GenBank files (*.gbk or *.gbff)

It will parse the organism name, and add the first two words of the organism name to the filename.

Note- expects folders inside of folders, like NCBI genomes. Needs to be
taught to work on a single folder

Also needs to be trained to remove periods from sp. to avoid sp..gbk
"""

import os
from pathlib import Path
from Bio import SeqIO

# Set the parent directory containing the folders
parent_dir = "/Users/daniel/Projects/wirehunter/all_3354_delta_wirehunter_results/folder/"  # <-- Update this path

# Loop through each folder in the parent directory
for folder in Path(parent_dir).iterdir():
    
    if not folder.is_dir():
        continue

    # Find all GenBank files (*.gbk or *.gbff) in the folder
    for gbk_file in folder.glob("*.gb*"):
        # Parse the GenBank file to get the organism name
        with open(gbk_file) as handle:
            record = next(SeqIO.parse(handle, "genbank"))
            organism = record.annotations.get("organism", "")
            print(f"Found {organism}")
        # Extract the first two words of the organism name
        organism_words = organism.split()
        if len(organism_words) >= 2:
            organism_prefix = f"{organism_words[0]}_{organism_words[1]}"
            print(f"Renaming {organism_words[0]}_{organism_words[1]}")
        else:
            organism_prefix = organism.replace(" ", "_")  # fallback

        # Construct the new filename
        new_filename = f"{gbk_file.stem}_{organism_prefix}{gbk_file.suffix}"
        new_path = gbk_file.with_name(new_filename)

        # Rename the file
        gbk_file.rename(new_path)
        print(f"✅ Renamed {gbk_file.name} → {new_filename}")

print("✅ All done!")