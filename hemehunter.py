#!/usr/bin/env python3

import argparse
import os
import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

#
# hemehunter - read gbk files, verify they have translated regions, and scan for CXXCH motifs
#  In proteins with more than 3 motifs, scan for noncannonical CxCH, CxxxCH, and others
#   
#   Output TSV file and optional fasta than can be used for further annotation r
#   such as adding those features to a genbank file as misc_features with coordinates
#
#   by ChattyGPT 4o and Daniel Bond
#   
#   updated 5/16/2024
#
#
# Simulate command-line arguments when running in an environment like Spyder.
# comment out for command line operation

sys.argv = [
    "hemehunter.py",
#    "-f",
    "GSU_DSMZ_3_23_translations.gbk",
]

# #
#
####################################################################

def parse_arguments():
    parser = argparse.ArgumentParser(
    description=(
        "Identify multiheme c-type cytochromes (≥3 CXXCH motifs, including noncannonical) from a GenBank file.\n\n"
        "This script scans translated CDS features and selects those with at least 3 canonical CXXCH motifs.\n"
        "It calculates motif counts, molecular weight, and motif density, and outputs a summary TSV.\n"
        "Use the -f flag to also save a FASTA file of the cytochrome amino acid sequences.\n\n"
        "It will automatically translate any untranslated regions, including truncated genes. If you don't\n"
        "want to annotate any new regions, just ignore the new file created ending in _translated."
    ),
    formatter_class=argparse.RawTextHelpFormatter
)
    parser.add_argument("genbank_file", help="Input GenBank file (.gbk or .gbff)")
    parser.add_argument(
        "-o", "--output", help="Output TSV file (default: <input>_cytochromes.tsv)"
    )
    parser.add_argument(
        "-f", "--fasta", action="store_true", help="Also output a FASTA file of the cytochrome protein sequences."
    )
    return parser.parse_args()


def has_translation(feature):
    """Returns True if feature has a valid /translation qualifier."""
    if "translation" not in feature.qualifiers:
        return False
    val = feature.qualifiers.get("translation", [""])[0]
    # print("🔍 Raw translation:", repr(val[:60]), "...")  # print first 60 characters
    return (
        isinstance(feature.qualifiers["translation"], list)
        and len(feature.qualifiers["translation"]) > 0
        and val.strip() != ""
    )

def ensure_translations(genbank_file):
    """
    Ensures all CDS features in the GenBank file have /translation qualifiers.
    If any are missing, translate and save new GenBank file with _translations.gbk suffix.
    Returns path to the usable GenBank file.
    """
    updated_records = []
    needs_update = False

    for record in SeqIO.parse(genbank_file, "genbank"):
        new_features = []
        for feature in record.features:
            if feature.type == "CDS":
                locus = feature.qualifiers.get("locus_tag", ["?"])[0]
                # print(f"\n[CHECK] locus_tag: {locus}")

                if has_translation(feature):
                    pass
                    #print(f"[OK] Translation already present for {locus}")
                else:
                    # print(f"[ADD] No valid translation found for {locus}, translating...")
                    try:
                        dna_seq = feature.extract(record.seq)
                        transl_table = int(feature.qualifiers.get("transl_table", [11])[0])
                        aa_seq = dna_seq.translate(table=transl_table, to_stop=True)
                        feature.qualifiers["translation"] = [str(aa_seq)]
                        needs_update = True  # ✅ Set only when we actually add one
                        # print(f"[SUCCESS] Translation added for {locus}")
                    except Exception as e:
                        print(f"[ERROR] Could not translate CDS at {locus}: {e}")

            new_features.append(feature)

        record.features = new_features
        updated_records.append(record)

    if needs_update:
        base = os.path.splitext(os.path.basename(genbank_file))[0]
        updated_file = base + "_translations.gbk"
        with open(updated_file, "w") as out_handle:
            SeqIO.write(updated_records, out_handle, "genbank")
        print(f"\n[INFO] CDS without translations found. New GenBank used in analysis written to: {updated_file}")
        return updated_file
    else:
        print("\n[INFO] All CDS features already have valid translations. No new file written.")
        return genbank_file

def check_genbank_file(genbank_file):
    """Validate input file contains CDS features with translatable sequences."""
    cds_found = False
    valid_translations = 0

    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                cds_found = True
                if has_translation(feature):
                    valid_translations += 1
                else:
                    try:
                        dna_seq = feature.extract(record.seq)
                        aa_seq = dna_seq.translate(table=1, to_stop=True)
                        if aa_seq:
                            valid_translations += 1
                    except Exception:
                        pass

    if not cds_found:
        raise ValueError("No CDS features found in GenBank file.")

    if valid_translations == 0:
        raise ValueError("No CDS features contain translatable sequences.")

    return True

####################################################################

def find_cxxch_motifs(protein_seq):
    """Return a dict with all CXXCH-like motif counts."""
    motifs = {
        "CXXCH": re.findall(r"C..CH", protein_seq),
        "CXXXCH": re.findall(r"C...CH", protein_seq),
        "CXCH": re.findall(r"C.CH", protein_seq),
        "CX14CH": re.findall(r"C.{10,14}CH", protein_seq)
    }
    return motifs

def extract_cytochrome_features(genbank_file):
    """
    Extracts CDS features with ≥3 CXXCH motifs and returns a list of dicts.
    Includes mapped start/end of first/last motif.
    """
    cytochromes = []

    motif_patterns = {
        "CXXCH": r"C..CH",
        "CXXXCH": r"C...CH",
        "CXCH": r"C.CH",
        "CX14CH": r"C.{10,14}CH"
    }

    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type != "CDS" or not has_translation(feature):
                continue

            aa_seq = feature.qualifiers["translation"][0]

            motif_locs = {name: [m.start() for m in re.finditer(pattern, aa_seq)]
                          for name, pattern in motif_patterns.items()}

            if len(motif_locs["CXXCH"]) < 3:
                continue  # Skip if fewer than 3 canonical motifs

            motif_positions = []
            for positions in motif_locs.values():
                motif_positions.extend(positions)

            if not motif_positions:
                continue

            motif_positions.sort()
            first_aa = motif_positions[0]
            last_aa = motif_positions[-1]

            cds_start = int(feature.location.start)
            cds_end = int(feature.location.end)
            strand = feature.location.strand

            if strand == 1:
                motif_start_genomic = cds_start + 3 * first_aa
                motif_end_genomic = cds_start + 3 * last_aa + 5
            else:
                motif_start_genomic = cds_end - 3 * last_aa - 5
                motif_end_genomic = cds_end - 3 * first_aa

            motif_count = sum(len(v) for v in motif_locs.values())
            odd_motifs = motif_count - len(motif_locs["CXXCH"])

            try:
                mw_kDa = molecular_weight(aa_seq, seq_type="protein") / 1000.0
            except Exception:
                mw_kDa = 0.0

            motif_density = motif_count / mw_kDa if mw_kDa > 0 else 0.0

            cytochromes.append({
                "locus_tag": feature.qualifiers.get("locus_tag", ["?"])[0],
                "product": feature.qualifiers.get("product", ["?"])[0],
                "gene": feature.qualifiers.get("gene", [""])[0],
                "motifs_total": motif_count,
                "motifs_noncanonical": odd_motifs,
                "kDa": round(mw_kDa, 2),
                "motifs_per_kDa": round(motif_density, 3),
                "protein_id": feature.qualifiers.get("protein_id", ["?"])[0],
                "start": min(motif_start_genomic, motif_end_genomic),
                "end": max(motif_start_genomic, motif_end_genomic),
                "aa_sequence": aa_seq,
            })

    return cytochromes

def write_cytochrome_tsv(features, output_path):
    """
    Writes tab-delimited tsv-style summary of multiheme cytochromes.
    """
    with open(output_path, "w") as f:
       header = [
    "locus_tag", "product", "gene", "motifs_total", "motifs_noncanonical",
    "kDa", "motifs_per_kDa", "protein_id", "start", "end"
       ]
       f.write("\t".join(header) + "\n")
       for feat in features:
            row = [str(feat.get(col, "")) for col in header]
            f.write("\t".join(row) + "\n")
    print(f"[INFO] Wrote {len(features)} cytochrome annotations to {output_path}")

def write_fasta(features, output_path):
    """
    Writes amino acid sequences of identified cytochromes to a FASTA file.
    """
    with open(output_path, "w") as f:
        for feat in features:
            header = f">{feat['locus_tag']}|{feat['product']}|{feat['protein_id']}|{feat['motifs_total']} hemes"
            seq = feat["aa_sequence"]
            f.write(header + "\n")
            for i in range(0, len(seq), 70):  # Wrap lines at 70 characters
                f.write(seq[i:i+70] + "\n")
    print(f"[INFO] Wrote FASTA file to {output_path}")
    
####################################################################    

if __name__ == "__main__":
    args = parse_arguments()

    # Determine output file name
    if args.output:
        output_file = args.output
    else:
        base = os.path.splitext(os.path.basename(args.genbank_file))[0]
        output_file = base + "_cytochromes.tsv"

    try:
        check_genbank_file(args.genbank_file)
        usable_file = ensure_translations(args.genbank_file)
        print(f"[INFO] GenBank file ready: {usable_file}")
    except ValueError as e:
        print(f"[ERROR] {e}")
        exit(1)
        
    # Process GenBank and extract cytochrome features

cyto_features = extract_cytochrome_features(usable_file)
write_cytochrome_tsv(cyto_features, output_file)# Continue processing from usable_file in later modules

if args.fasta:
    fasta_name = os.path.splitext(output_file)[0] + ".fasta"
    write_fasta(cyto_features, fasta_name)