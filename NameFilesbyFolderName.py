#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 this looks at a folder *full of folders* named by their Acession Number, like GCF_0238273.1

   inside each folder are files with extensions.
   
   This will just remove the filename and

   replace it with the folder's name, plus the original acession. So if the files are all

  genome.gbk and genome.faa, they become GCF_0238273.1.gbk and GCF_0238273.1.faa
  
  If you just have a gbff file, make sure it is in a folder

"""

import os




# There is no argparsing, just run from Spyder

parent_dir = "/Users/daniel/Projects/wirehunter/all_3354_delta_wirehunter_results/4_betasearch/"  

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