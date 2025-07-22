#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
when you have a folder of gbk or gbff files

you need to extract all the locus and organism names

so you can link GCF or acession numbers to contigs

"""


import os
import re
import csv
from pathlib import Path
from Bio import SeqIO

# Input folder to scan
input_folder = "/Users/daniel/Projects/Downloaded_genomes/Thermodesulfobacteriota_Deltas_RefSeq_848/Delta_gbff/"
output_csv = "/Users/daniel/Projects/Downloaded_genomes/Thermodesulfobacteriota_Deltas_RefSeq_848/Delta_gbff/genbank_metadata.csv"

# Pattern to extract accession (e.g., GCF_010646885.2)
accession_pattern = re.compile(r'(GCF_\d+\.\d+)')

# Collect metadata
records = []

for root, _, files in os.walk(input_folder):
    for file in files:
        if file.endswith(('.gbk', '.gbff')):
            match = accession_pattern.search(file)
            if not match:
                print(f"⚠️ Could not extract accession from filename: {file}")
                continue
            accession = match.group(1)
            filepath = os.path.join(root, file)
            try:
                with open(filepath, "r") as handle:
                    # Take the first record (assumes GenBank file has at least one record)
                    record = next(SeqIO.parse(handle, "genbank"))

                    locus = record.name  # LOCUS line
                    definition = record.description  # DEFINITION line
                    organism = record.annotations.get("organism", "Unknown")

                    records.append([accession, locus, definition, organism])
            except Exception as e:
                print(f"❌ Error parsing {file}: {e}")

# Write to CSV
with open(output_csv, "w", newline='') as out:
    writer = csv.writer(out)
    writer.writerow(["Accession", "Locus", "Definition", "Organism"])
    writer.writerows(records)

print(f"✅ Metadata written to: {output_csv}")