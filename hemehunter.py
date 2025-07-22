#!/usr/bin/env python3

import argparse
import os
import sys
import re
from Bio import SeqIO
from Bio.SeqUtils import molecular_weight
import warnings
from Bio import BiopythonWarning

# # #
#  Simulate command-line arguments when running in an environment like Spyder.
#  Comment out for command line operation

# sys.argv = [
#     "hemehunter.py",
#     "/Users/daniel/Desktop/ReturnoftheChrome/Desulfuro_results/Desulfo_Signal/gbks/",
#     "-f",
#     "-m 2",
#     "--force"
# ]

####################################################################
# Parser logic
def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "****************************************************************************************\n"
            "HemeHunter 3000 identifies multiheme c-type cytochromes in a GenBank file or directory.\n\n"
            "This script scans translated CDS features and selects those with (default) at least 3,\n"
            "or the threshold set with -m, to find that many CXXCH motifs before looking for others.\n"
            "The noncannonical list is CXCH, CXXXCH, and C (10,14) CH, and avoids overcounting errors.\n"
            "It calculates motif counts, molecular weight, and motif density, and outputs a summary TSV.\n"
            "Use the -f flag to also save a FASTA file of the cytochrome amino acid sequences.\n\n"
            "If you need translations, use -t to auto-translate missing regions and write _translations.gbk.\n"
            "It will autodetect directories or single files for batch processing.\n"
            "****************************************************************************************\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("genbank_file", help="Input GenBank file (.gbk or .gbff) or a folder for batch mode")
    parser.add_argument("-o", "--output", help="Output TSV file (default: <input>_cytochromes.tsv)")
    parser.add_argument("-f", "--fasta", action="store_true", help="Optional: output a FASTA file of cytochrome proteins.")
    parser.add_argument("-m", "--min_motifs", type=int, default=3, help="Minimum number of CXXCH motifs (default: 3)")
    parser.add_argument("-t", "--translations", action="store_true", help="If set, translate missing CDS features.")
    parser.add_argument("--force", action="store_true", help="Overwrite existing output files if present.")
    return parser

def parse_arguments():
    parser = build_parser()
    return parser.parse_args()

####################################################################
# Other functions unchanged
def has_translation(feature):
    if "translation" not in feature.qualifiers:
        return False
    val = feature.qualifiers.get("translation", [""])[0]
    return isinstance(feature.qualifiers["translation"], list) and len(feature.qualifiers["translation"]) > 0 and val.strip() != ""

def ensure_translations(genbank_file, overwrite=False):
    updated_records = []
    needs_update = False

    for record in SeqIO.parse(genbank_file, "genbank"):
        new_features = []
        for feature in record.features:
            if feature.type == "CDS":
                if not has_translation(feature):
                    try:
                        dna_seq = feature.extract(record.seq)
                        transl_table = int(feature.qualifiers.get("transl_table", [11])[0])
                        trimmed_len = len(dna_seq) - (len(dna_seq) % 3)
                        dna_seq = dna_seq[:trimmed_len]
                        aa_seq = dna_seq.translate(table=transl_table, to_stop=True)
                        feature.qualifiers["translation"] = [str(aa_seq)]
                        needs_update = True
                    except Exception as e:
                        print(f"[ERROR] Could not translate CDS: {e}")
            new_features.append(feature)
        record.features = new_features
        updated_records.append(record)

    if needs_update:
        base = os.path.splitext(os.path.basename(genbank_file))[0]
        input_dir = os.path.dirname(os.path.abspath(genbank_file))
        updated_file = os.path.join(input_dir, base + "_translations.gbk")

        if os.path.exists(updated_file) and not overwrite:
            print(f"[WARNING] {updated_file} exists. Use --force to overwrite. Skipping translation.")
            return genbank_file

        with open(updated_file, "w") as out_handle:
            SeqIO.write(updated_records, out_handle, "genbank")
        print(f"[INFO] Translations added. New file: {updated_file}")
        return updated_file
    else:
        print("[INFO] All CDS features have valid translations.")
        return genbank_file

def check_genbank_file(genbank_file):
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
        raise ValueError("No CDS features found.")
    if valid_translations == 0:
        raise ValueError("No CDS features with translatable sequences.")
    return True

def extract_cytochrome_features(genbank_file, min_motifs=3):
    cytochromes = []

    noncanonical_patterns = {
        "CXXXCH": re.compile(r'(C[^CH]{3})(CH)'),
        "CXCH": re.compile(r'(C[^CH]{1})(CH)'),
        "CX14CH": re.compile(r'(C.{10,14}?)(CH)')
    }

    def greedy_find_cxxch(seq):
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

            canon_matches, used_ch = greedy_find_cxxch(aa_seq)
            motif_positions = [start for start, _ in canon_matches]
            noncanon_counts = {"CXCH": 0, "CXXXCH": 0, "CX14CH": 0}

            for name, pattern in noncanonical_patterns.items():
                for m in pattern.finditer(aa_seq):
                    h_index = m.start(2)
                    ch_pos = (h_index, h_index + 1)
                    if ch_pos not in used_ch:
                        motif_positions.append(m.start())
                        used_ch.add(ch_pos)
                        noncanon_counts[name] += 1

            if len(canon_matches) + sum(noncanon_counts.values()) < min_motifs:
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

            try:
                mw_kDa = molecular_weight(aa_seq, seq_type="protein") / 1000.0
            except Exception:
                mw_kDa = 0.0

            motif_density = (len(canon_matches) + sum(noncanon_counts.values())) / mw_kDa if mw_kDa > 0 else 0.0

            cytochromes.append({
                "locus_tag": feature.qualifiers.get("locus_tag", ["?"])[0],
                "product": feature.qualifiers.get("product", ["?"])[0],
                "gene": feature.qualifiers.get("gene", [""])[0],
                "motifs_total": len(canon_matches) + sum(noncanon_counts.values()),
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

def write_cytochrome_tsv(features, output_path, overwrite=False):
    if os.path.exists(output_path) and not overwrite:
        print(f"[WARNING] {output_path} exists. Use --force to overwrite.")
        return
    if not features:
        print(f"[INFO] No cytochromes found for {output_path}")
        return
    with open(output_path, "w", newline="") as f:
        header = [
            "locus_tag", "product", "gene", "motifs_total",
            "motifs_CXCH", "motifs_CXXXCH", "motifs_CX10_14CH",
            "kDa", "motifs_per_kDa", "protein_id", "start", "end"
        ]
        f.write("\t".join(header) + "\n")
        for feat in features:
            row = [str(feat.get(col, "")) for col in header]
            f.write("\t".join(row) + "\n")
    print(f"[INFO] Wrote {len(features)} annotations to {output_path}")

def write_fasta(features, output_path):
    with open(output_path, "w") as f:
        for feat in features:
            header = f">{feat['locus_tag']}|{feat['product']}|{feat['protein_id']}|{feat['motifs_total']} hemes"
            seq = feat["aa_sequence"]
            f.write(header + "\n")
            for i in range(0, len(seq), 70):
                f.write(seq[i:i+70] + "\n")
    print(f"[INFO] Wrote FASTA to {output_path}")

def run_on_one_genbank_file(gbk_path, args):
    try:
        check_genbank_file(gbk_path)
        usable_file = ensure_translations(gbk_path, overwrite=args.force) if args.translations else gbk_path
        input_path = os.path.abspath(usable_file)
        input_dir = os.path.dirname(input_path)
        base_name = os.path.splitext(os.path.basename(input_path))[0]

        output_file = args.output
        if not output_file:
            output_file = os.path.join(input_dir, f"{base_name}_cytochromes.tsv")
        elif os.path.isdir(output_file):
            output_file = os.path.join(output_file, f"{base_name}_cytochromes.tsv")

        fasta_name = os.path.join(input_dir, f"{base_name}_cytochromes.fasta")

        cyto_features = extract_cytochrome_features(usable_file, min_motifs=args.min_motifs)
        write_cytochrome_tsv(cyto_features, output_file, overwrite=args.force)
        if args.fasta:
            write_fasta(cyto_features, fasta_name)
    except ValueError as e:
        print(f"[ERROR] {gbk_path}: {e}")

####################################################################    
# Main entry point
def main():
    args = parse_arguments()
    input_path = os.path.abspath(args.genbank_file)
    if os.path.isdir(input_path):
        print(f"[INFO] Batch mode: scanning {input_path}...")
        for fname in os.listdir(input_path):
            if fname.endswith((".gb", ".gbk", ".gbff")):
                full_path = os.path.join(input_path, fname)
                print(f"📄 Processing: {fname}")
                run_on_one_genbank_file(full_path, args)
    else:
        run_on_one_genbank_file(input_path, args)

if __name__ == "__main__":
    main()