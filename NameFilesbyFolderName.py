#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 18 21:58:34 2025

@author: daniel
"""

import os


# this looks at a folder full of folders named by their Acession Number, like GCF_0238273.1
#   inside each folder are files with extensions, this will just remove the filename and
#   replace it with the folder's name, plus the original acession. So if the files are all
#   genome.gbk and genome.faa, they become GCF_0238273.1.gbk and GCF_0238273.1.faa

# There is no argparsing, just run from Spyder

parent_dir = "/Users/daniel/Desktop/ReturnoftheChrome/DesulfuromonadalesRefSeq"  

# replace with the path to your folder containing the folders

for accession in os.listdir(parent_dir):
    folder_path = os.path.join(parent_dir, accession)
    if not os.path.isdir(folder_path):
        continue

    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)

        # Extract original extension
        ext = os.path.splitext(filename)[1]

        # Create new filename
        new_name = f"{accession}{ext}"
        new_path = os.path.join(folder_path, new_name)

        # Rename file
        os.rename(file_path, new_path)
        print(f"✅ Renamed {filename} → {new_name}")