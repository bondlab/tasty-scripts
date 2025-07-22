#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
when you have a list of GCF or GCA acessions

you need to go to NCBI and get their metadata

so you can link various loci and accesions to organisms

because you should have done this earlier

"""
from Bio import Entrez
import csv
import time
import re

Entrez.email = "dbond@umn.edu"  # Replace with your actual email

# Input file with one GCF or GCA accession per line
# Will populate a CSV with name and other NCBI data

input_file = "/Users/daniel/Projects/Downloaded_genomes/All_delta_GenBank_3354/gb_acessions_toget.txt"
output_file = "/Users/daniel/Projects/Downloaded_genomes/All_delta_GenBank_3354/gb_acessions_metadata.csv"

# Read accession numbers from file
with open(input_file) as f:
    accessions = [line.strip() for line in f if line.strip()]

# Open output CSV
with open(output_file, "w", newline='') as out_csv:
    writer = csv.writer(out_csv)
    writer.writerow(["Accession", "Organism", "Species", "TaxID", "Status", "Assembly", "Reference"])

    for acc in accessions:
        try:
            # Step 1: Search NCBI Assembly to get UID
            search = Entrez.esearch(db="assembly", term=acc)
            record = Entrez.read(search)
            search.close()

            if not record["IdList"]:
                print(f"⚠️ Accession {acc} not found.")
                continue

            uid = record["IdList"][0]

            # Step 2: Fetch the summary metadata
            summary = Entrez.esummary(db="assembly", id=uid, report="full")
            result = Entrez.read(summary)
            summary.close()

            doc = result['DocumentSummarySet']['DocumentSummary'][0]

            organism_raw = doc.get('Organism', '')
            organism = re.sub(r'\s*\([^)]*\)', '', organism_raw).strip() #take out (bacteria), etc
            species = doc.get('SpeciesName')
            assembly = doc.get('AssemblyName', '')
            taxid = doc.get('Taxid', '')
            status = doc.get('AssemblyStatus', '')
            ref = doc.get('RefSeq_category', '')

            writer.writerow([acc, organism, species, taxid, status, assembly, ref])
            print(f"✅ Retrieved {acc}: {organism}")

            time.sleep(0.4)  # NCBI rate limit

        except Exception as e:
            print(f"❌ Error fetching {acc}: {e}")