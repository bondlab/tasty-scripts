#!/usr/bin/env python3

import os
import sys
from pathlib import Path


sys.argv = [
    "extract_gbk_info.py",
    "/Users/daniel/Projects/wirehunter/all_3354_delta_wirehunter_results/unique_gbk_with_wires/",
    "/Users/daniel/Projects/wirehunter/all_3354_delta_wirehunter_results/uniquque_gbk_metadata.tsv",
]



def extract_locus_and_organism(gbk_path):
    locus = None
    organism = None
    with open(gbk_path, 'r') as file:
        for line in file:
            if line.startswith("LOCUS"):
                locus = line.split()[1]
            elif line.strip().startswith("ORGANISM"):
                organism = line.strip().split(" ", 1)[1]
            if locus and organism:
                break
    return locus, organism

def main(input_dir, output_tsv):
    input_path = Path(input_dir)
    gbk_files = list(input_path.glob("*.gbk"))

    with open(output_tsv, 'w') as out:
        out.write("LOCUS\tORGANISM\tFILENAME\n")
        for gbk_file in gbk_files:
            locus, organism = extract_locus_and_organism(gbk_file)
            if locus and organism:
                name = gbk_file.stem  # removes '.gbk'
                out.write(f"{locus}\t{organism}\t{name}\n")
            else:
                print(f"⚠️ Could not extract data from {gbk_file}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python extract_gbk_info.py <input_folder> <output.tsv>")
    else:
        main(sys.argv[1], sys.argv[2])