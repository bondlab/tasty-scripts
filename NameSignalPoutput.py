#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Renames 'output.gff3' inside each accession-named folder to '{accession}_signaloutput.gff3'
"""

import os

# this looks at a folder full of folders named by their Acession Number, like GCF_0238273.1
#   inside each folder is an output.gff3 file. This will rename the output file so 
#   output.gff3 becomes GCF_0238273.1_signaloutput.gff3

# There is no argparsing, just run from Spyder

parent_dir = "/Users/daniel/Desktop/ReturnoftheChrome/SignalP"  # adjust as needed

for accession in os.listdir(parent_dir):
    folder_path = os.path.join(parent_dir, accession)
    if not os.path.isdir(folder_path):
        continue

    old_file = os.path.join(folder_path, "output.gff3")
    if os.path.isfile(old_file):
        new_file = os.path.join(folder_path, f"{accession}_signaloutput.gff3")
        os.rename(old_file, new_file)
        print(f"✅ Renamed output.gff3 → {accession}_signaloutput.gff3")
    else:
        print(f"⚠️ No output.gff3 found in {accession}")