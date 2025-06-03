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
import warnings
from Bio import BiopythonWarning

#
# hemehunter - read gbk files, verify they have translated regions, and scan for CXXCH motifs
#  In proteins with more than 3 motifs, scan for noncannonical CxCH, CxxxCH, and others
#   
#   Output TSV file and optional fasta than can be used for further annotation r
#   such as adding those features to a genbank file as misc_features with coordinates
#
#   by ChattyGPT 4o and Daniel Bond
#   
#   updated 5/26/2024

#
#
#
# Simulate command-line arguments when running in an environment like Spyder.
# comment out for command line operation

sys.argv = [
    "hemehunter.py",
    "//Users/daniel/Desktop/ReturnoftheChrome/GeobacteraceaeRefSeq/Geo_gbffs",
    "-f",
    "-m 2",
    "--force"  
]

# # #
# 

####################################################################

# Parse files and prepare translation if needed

def parse_arguments():
    parser = argparse.ArgumentParser(
    description=(
        "****************************************************************************************\n"
        "HemeHunter 3000 identifies multiheme c-type cytochromes in a GenBank file or directory of files.\n\n"
        "This script scans translated CDS features and selects those with (default) at least 3,\n"
        "or the threshold set with -m, to find that many CXXCH motifs before looking for others.\n"
        "The noncannonical list is CXCH, CXXXCH, and C (10,14) CH, and avoids overcounting errors.\n"
        "It calculates motif counts, molecular weight, and motif density, and outputs a summary TSV.\n"
        "Use the -f flag to also save a FASTA file of the cytochrome amino acid sequences.\n\n"
        "If you need translations, or want to look for missing genes, use the option -t \n"
        "which will automatically translate any untranslated regions, including truncated genes.\n"
        "It will then use the new _translated.gbk version of the file to call cytochromes\n"
        "\nIt will autodetect if you provide a directory or single file, and process appropriately\n"
        "****************************************************************************************\n"
    ),
    formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("genbank_file", help="Input GenBank file (.gbk or .gbff) or a folder for batch mode")
    parser.add_argument("-o", "--output", help="Output TSV file (default: <input>_cytochromes.tsv)")
    parser.add_argument("-f", "--fasta", action="store_true", help="Optional, output a FASTA file of the cytochrome protein sequences.")
    parser.add_argument("-m", "--min_motifs", type=int, default=3, help="Minimum number of CXXCH motifs to consider a protein a multiheme cytochrome (default: 3)")
    parser.add_argument("-t", "--translations", action="store_true", help="If set, translate CDS features missing /translation qualifiers and write a new _translations.gbk file.")
    parser.add_argument(
    "--force", action="store_true",
    help="If set, allows overwriting existing output files (_translations.gbk, .tsv, .fasta)."
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

def ensure_translations(genbank_file, overwrite=False):
    """
    Ensures all CDS features in the GenBank file have /translation qualifiers.
    If any are missing, translate and save new GenBank file with _translations.gbk suffix.
    Returns path to the usable GenBank file.
    """
    updated_records = []
    needs_update = False

    print(f"[DEBUG] overwrite = {overwrite}")

    for record in SeqIO.parse(genbank_file, "genbank"):
        new_features = []
        for feature in record.features:
            if feature.type == "CDS":
                locus = feature.qualifiers.get("locus_tag", ["?"])[0]
                if has_translation(feature):
                    pass
                else:
                    try:
                        dna_seq = feature.extract(record.seq)
                        transl_table = int(feature.qualifiers.get("transl_table", [11])[0])
                        trimmed_len = len(dna_seq) - (len(dna_seq) % 3)
                        dna_seq = dna_seq[:trimmed_len]
                        aa_seq = dna_seq.translate(table=transl_table, to_stop=True)
                        feature.qualifiers["translation"] = [str(aa_seq)]
                        needs_update = True
                    except Exception as e:
                        print(f"[ERROR] Could not translate CDS at {locus}: {e}")
            # ✅ Always add the feature (whether gene or CDS or other)
            new_features.append(feature)
    
        record.features = new_features
        updated_records.append(record)

    if needs_update:
        base = os.path.splitext(os.path.basename(genbank_file))[0]
        input_dir = os.path.dirname(os.path.abspath(genbank_file))
        updated_file = os.path.join(input_dir, base + "_translations.gbk")
        
        if os.path.exists(updated_file) and not overwrite:
            print(f"[WARNING] {updated_file} already exists. Use --force to overwrite. Skipping translation step.")
            return genbank_file
    
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

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)

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

# Collect features. Can add other noncannonical motifs if desired

def extract_cytochrome_features(genbank_file, min_motifs=3):
    """
    Extracts CDS features with ≥min_motifs CXXCH motifs or total motifs (non-overlapping).
    Tracks motif types and genomic location boundaries.
    
    Once a CH is counted as part of one motif, it will not be allowed to be part of another, to
    avoid accidential overcounting errors when motifs are close together, as in C 10-14 CH,
    or have extra C's, as in CCXCH, etc, which could be be double-counted as CXXCH and CXCH'
    """

    cytochromes = []

    noncanonical_patterns = {
        "CXXXCH": re.compile(r'(C[^CH]{3})(CH)'),
        "CXCH": re.compile(r'(C[^CH]{1})(CH)'),
        "CX14CH": re.compile(r'(C.{10,14}?)(CH)')  # use non-greedy match
    }
    
    def greedy_find_cxxch(seq):
        """Greedy left-to-right non-overlapping CXXCH finder."""
        pattern = re.compile(r"C..CH")
        used_ch = set()
        matches = []
        pos = 0

        while pos <= len(seq) - 5:
            m = pattern.search(seq, pos)
            if not m:
                break
            ch_pos = (m.start() + 3, m.start() + 4)
            if ch_pos not in used_ch:
                matches.append((m.start(), m.end()))
                used_ch.add(ch_pos)
                pos = m.end()
            else:
                pos = m.start() + 1

        return matches, used_ch

    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type != "CDS" or "translation" not in feature.qualifiers:
                continue

            aa_seq = feature.qualifiers["translation"][0]
            motif_positions = []
            # odd_motif_count = 0  # old way of counting motifs


            # === Stage 1: Greedy canonical CXXCH search === 
            canon_matches, used_ch = greedy_find_cxxch(aa_seq)
            motif_positions = [start for start, _ in canon_matches]
            noncanon_counts = {"CXCH": 0, "CXXXCH": 0, "CX14CH": 0}
            
           # === Stage 2: Always scan for noncanonical matches, skipping overlapping CH ===
            for name, pattern in noncanonical_patterns.items():
                for m in pattern.finditer(aa_seq):
                    h_index = m.start(2)  # Start of the "CH" match
                    ch_pos = (h_index, h_index + 1)
            
                    if ch_pos not in used_ch:
                        motif_positions.append(m.start())
                        used_ch.add(ch_pos)
                        noncanon_counts[name] += 1
            
            # === Threshold filter AFTER all motif types counted ===
            if len(canon_matches) + sum(noncanon_counts.values()) < min_motifs:
                continue  # still below threshold; skip this CDS
                
            # === Finalize total motif count and boundaries ===
            total_motifs = len(canon_matches) + sum(noncanon_counts.values())
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
            
            # === Compute MW and density ===
            try:
                mw_kDa = molecular_weight(aa_seq, seq_type="protein") / 1000.0
            except Exception:
                mw_kDa = 0.0
            
            motif_density = total_motifs / mw_kDa if mw_kDa > 0 else 0.0
            
            # === Store result ===
            cytochromes.append({
                "locus_tag": feature.qualifiers.get("locus_tag", ["?"])[0],
                "product": feature.qualifiers.get("product", ["?"])[0],
                "gene": feature.qualifiers.get("gene", [""])[0],
                "motifs_total": total_motifs,
                "motifs_CXCH": noncanon_counts["CXCH"],
                "motifs_CXXXCH": noncanon_counts["CXXXCH"],
                "motifs_CX10_14CH": noncanon_counts["CX14CH"],
                "kDa": round(mw_kDa, 2),
                "motifs_per_kDa": round(motif_density, 3),
                "protein_id": feature.qualifiers.get("protein_id", ["?"])[0],
                "start": min(motif_start_genomic, motif_end_genomic),
                "end": max(motif_start_genomic, motif_end_genomic),
                "aa_sequence": aa_seq,
            })

    return cytochromes

####################################################################

# Write files of collected features

def write_cytochrome_tsv(features, output_path, overwrite=False):
    if os.path.exists(output_path) and not overwrite:
        print(f"[WARNING] {output_path} already exists. Use --force to overwrite.")
        return

    if not features:
        print(f"[INFO] No cytochromes found. Skipping TSV write for {output_path}")
        return

    print(f"[DEBUG] Writing TSV with {len(features)} features to {output_path}")

    with open(output_path, "w") as f:
       header = [
    "locus_tag", "product", "gene", "motifs_total",
    "motifs_CXCH", "motifs_CXXXCH", "motifs_CX10_14CH",
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

# Main loop that processes one genbank file at a time

def run_on_one_genbank_file(gbk_path, args):
    try:
        check_genbank_file(gbk_path)

        if args.translations:
            usable_file = ensure_translations(gbk_path, overwrite=args.force)
        else:
            usable_file = gbk_path

        print(f"[INFO] GenBank file ready: {usable_file}")

        # Now define output paths based on the final usable_file
        input_path = os.path.abspath(usable_file)
        input_dir = os.path.dirname(input_path)
        base_name = os.path.splitext(os.path.basename(input_path))[0]

        # Construct default output paths
        if args.output:
            if os.path.isdir(args.output):
                output_file = os.path.join(args.output, f"{base_name}_cytochromes.tsv")
            else:
                output_file = args.output
        else:
            output_file = os.path.join(input_dir, f"{base_name}_cytochromes.tsv")

        fasta_name = os.path.join(input_dir, f"{base_name}_cytochromes.fasta")

        # Extract and write results
        cyto_features = extract_cytochrome_features(usable_file, min_motifs=args.min_motifs)
        write_cytochrome_tsv(cyto_features, output_file, overwrite=args.force)

        if args.fasta:
            write_fasta(cyto_features, fasta_name)

    except ValueError as e:
        print(f"[ERROR] {gbk_path}: {e}")

####################################################################    

# main is primarily deciding if this is batch or single mode, and sending files to run_on_one_genbank_file

if __name__ == "__main__":
    args = parse_arguments()
    input_path = os.path.abspath(args.genbank_file)

    if os.path.isdir(input_path):
        print(f"[INFO] Batch mode: scanning {input_path} for GenBank files...")
        for fname in os.listdir(input_path):
            if fname.endswith((".gb", ".gbk", ".gbff")):
                full_path = os.path.join(input_path, fname)
                print(f"📄 Processing file: {fname}")
                run_on_one_genbank_file(full_path, args)
    else:
        run_on_one_genbank_file(input_path, args)