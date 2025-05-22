#!/usr/bin/env python3

#
#   fastfasta - extract amino acid sequences from a genbank file

import argparse
import os
from Bio import SeqIO

def extract_proteins_to_fasta(gbk_file):
    base = os.path.splitext(gbk_file)[0]
    output_file = base + ".fasta"

    count = 0
    with open(output_file, "w") as out_f:
        for record in SeqIO.parse(gbk_file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS" and "translation" in feature.qualifiers:
                    locus = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                    product = feature.qualifiers.get("product", ["hypothetical protein"])[0]
                    aa_seq = feature.qualifiers["translation"][0]
                    header = f">{locus} {product}"
                    out_f.write(f"{header}\n{aa_seq}\n")
                    count += 1
    print(f"✅ Wrote {count} protein sequences to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Extract protein sequences from a GenBank file and write to a FASTA file.")
    parser.add_argument("genbank_file", help="Input GenBank file (.gb or .gbk)")
    args = parser.parse_args()
    extract_proteins_to_fasta(args.genbank_file)

if __name__ == "__main__":
    main()