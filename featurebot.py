#!/usr/bin/env python3

import argparse
import os
import re
import sys
import csv
import glob
from pathlib import Path
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

#
# Featurebot - read GFF files resulting from various prediction tools
#   and add those features to a genbank file as misc_features with coordinates
#
#   by ChattyGPT 4o and Daniel Bond
#   
#   updated 5/23/2024
#
#
#   to do: clean up the expectation that the input is a -gff and change all 
#   occurences to -file or something. won't change behavior but reflects types of files
#
#
#
# Simulate command-line arguments when running in an environment like Spyder.
# comment out for command line operation


sys.argv = [
    "featurebot.py",
    "-gbkdir", "/Users/daniel/Desktop/ReturnoftheChrome/Results to map onto Geo/Heme_Geo_analysis/gbks_to_map",
    "-featuredir", "/Users/daniel/Desktop/ReturnoftheChrome/Results to map onto Geo/Heme_Geo_analysis/HemeResults_tomap",
    "-multiheme",
    "--force"
]

# 
#
####################################################################

def parse_arguments():
    """
    Parse command-line arguments for GenBank and GFF or TSV files,
    including analysis type and automatic output filename generation.
    """
    parser = argparse.ArgumentParser(
        description="Update a GenBank file with output from signal peptide, beta-barrel, or multiheme analyses.\n"
        "Can run in single mode (use -gff and -gbk for feature and genome file), or batch (specify folders).\n"
        "\n"
        "In batch mode, will match features to genomes if acession number is filename prefix, like GCF_000816575.1_regionoutput.gff3."
    )
    
    parser.add_argument(
        "-gff", "--gff_file",
        required=False,
        help="Path to the input file containing predicted features. Accepts gff3 output from SignalP 6.0, DTU DeepTHMM, and tsv from hemehunter"
    )
    
    parser.add_argument(
        "-gbk", "--genbank_file",
        required=False,
        help="Path to the input GenBank (.gb or .gbk or .gbff) file."
    )
    
    parser.add_argument(
        "-gbkdir", "--genbank_directory",
        required=False,
        help="Directory containing GenBank (.gbk, .gb, or .gbff) files for batch processing."
    )
    
    parser.add_argument(
        "-featuredir", "--feature_directory",
        required=False,
        help="Directory containing GFF/TSV feature files with matching prefixes. Default is to match locus tags, use -proteinid to use protein_id."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force overwrite of existing output files in batch mode."
    )
    
    parser.add_argument(
    "-proteinid",
    action="store_true",
    help="Use protein_id instead of locus_tag to match CDS features to annotations in -signal or -beta mode."
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-signal", action="store_true", help="Parse signal peptide predictions from SignalP 6")
    group.add_argument("-beta", action="store_true", help="Parse beta-barrel predictions from DeepTHMM")
    group.add_argument("-multiheme", action="store_true", help="Parse multiheme cytochrome predictions")
    
    parser.add_argument(
        "-o", "--output",
        required=False,
        default=None,
        help="Output GenBank file with updated annotations (default: original name + _signal.gbk, _beta.gbk, or _multiheme.gbk)"
    )

    args = parser.parse_args()
    
    # Only generate output name in single-file mode
    if args.output is None and args.genbank_file and args.gff_file:
        base, _ = os.path.splitext(os.path.basename(args.genbank_file))
        suffix = "_signal.gbk" if args.signal else "_beta.gbk" if args.beta else "_multiheme.gbk"
        args.output = base + suffix
    
    return args

####################################################################

def check_input_files(gff_path, gbk_path, args):
    """
    Check that the input files exist and have valid extensions, based on analysis type.
    """
    if not os.path.isfile(gff_path):
        sys.exit(f"❌ Error: file not found: {gff_path}")

    if args.multiheme:
        if not gff_path.lower().endswith(".tsv"):
            sys.exit(f"❌ Error: expected a .tsv file for multiheme analysis, got: {gff_path}")
    else:
        if not os.path.isfile(gbk_path) or not gbk_path.lower().endswith((".gb", ".gbk", ".gbff")):
            sys.exit(f"❌ Error: expected a .gff or .gff3 file for signal or beta analysis, got: {gff_path}")

    if not os.path.isfile(gbk_path) or not gbk_path.lower().endswith((".gb", ".gbk")):
        sys.exit(f"❌ Error: {gbk_path} is not a valid .gb or .gbk file.")

# Signal P parsing ###################################################################

def parse_signalp_gff(gff_file):
    """
    Parse a SignalP-style GFF file with non-standard format.
    Returns a dictionary keyed by locus_tag with a list of feature dicts.
    """
    annotations = {}

    with open(gff_file, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue  # skip comments and blank lines
            fields = line.strip().split("\t")
            if len(fields) < 5:
                continue  # skip malformed lines

            # First field: e.g., "GSU0015 protein description"
            id_and_desc = fields[0]
            locus_tag = id_and_desc.split()[0]
            description = " ".join(id_and_desc.split()[1:])
            source = fields[1]
            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])

            feature = {
                "type": feature_type,
                "start": start,
                "end": end,
                "source": source,
                "description": description
            }

            if locus_tag not in annotations:
                annotations[locus_tag] = []
            annotations[locus_tag].append(feature)

    return annotations

####################################################################

def extract_signal_peptides(annotation_dict):
    """
    For each locus tag, extract the start of the n-region and end of the h-region.
    Returns a dict: locus_tag → (n_start, h_end)
    """
    signal_peptides = {}
    for locus_tag, features in annotation_dict.items():
        n_start = None
        h_end = None
        for f in features:
            if f["type"] == "n-region":
                n_start = f["start"]
            elif f["type"] == "h-region":
                h_end = f["end"]
        if n_start and h_end:
            signal_peptides[locus_tag] = (n_start, h_end)
    return signal_peptides

####################################################################

def extract_lipid_cysteines(annotation_dict):
    """
    For each locus tag, extract the position of lipid-modified cysteine (aa coordinate).
    Returns a dict: locus_tag → residue_index
    """
    lipid_sites = {}
    for locus_tag, features in annotation_dict.items():
        for f in features:
            if f["type"] == "lipid-modified cysteine":
                if f["start"] == f["end"]:
                    lipid_sites[locus_tag] = f["start"]
    return lipid_sites

####################################################################

def extract_tat_motifs(annotation_dict):
    """
    For each locus tag, extract the start and end of the twin-arginine motif.
    Returns a dict: locus_tag → (start_aa, end_aa)
    """
    tat_motifs = {}
    for locus_tag, features in annotation_dict.items():
        for f in features:
            if f["type"] == "twin-arginine motif":
                tat_motifs[locus_tag] = (f["start"], f["end"])
    return tat_motifs

# TMHMM parsing ################################################################### 

def parse_beta_gff(gff_path):
    """
    Parse a beta-barrel topology GFF-style file and extract feature regions by locus tag.

    Returns:
        beta_features: dict mapping locus_tag -> list of dicts with 'type', 'start', and 'end'
    """
    beta_features = {}
    current_tag = None

    with open(gff_path, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("//"):
                continue  # Skip comments and separators

            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue  # skip lines that don't have enough fields

            locus_tag = parts[0]
            region_type = parts[1].strip()
            try:
                start = int(parts[2])
                end = int(parts[3])
            except ValueError:
                continue  # skip malformed lines
            
            locus_tag = parts[0]
            region_type = parts[1].strip()
            try:
                start = int(parts[2])
                end = int(parts[3])
            except ValueError:
                continue  # skip malformed lines
            locus_tag = parts[0]
            region_type = parts[1]
            try:
                start = int(parts[2])
                end = int(parts[3])
            except ValueError:
                continue  # Skip lines with non-integer coordinates

            if locus_tag not in beta_features:
                beta_features[locus_tag] = []

            beta_features[locus_tag].append({
                'type': region_type,
                'start': start,
                'end': end
            })

    return beta_features
  
# multiheme parsing ###################################################################

def parse_multiheme_tsv(gff_file):
    """
    Parses a multiheme cytochrome TSV file into a dictionary by locus_tag.

    Returns:
        A dict mapping locus_tag -> dict with product, gene, motifs, etc.
    """
    features = {}

    with open(gff_file, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            locus = row["locus_tag"]
            features[locus] = {
                "product": row["product"],
                "gene": row["gene"],
                "motifs_total": int(row["motifs_total"]),
                "motifs_CXCH": int(row["motifs_CXCH"]),
                "motifs_CXXXCH": int(row["motifs_CXXXCH"]),
                "motifs_CX10_14CH": int(row["motifs_CX10_14CH"]),
                "kDa": float(row["kDa"]),
                "motifs_per_kDa": float(row["motifs_per_kDa"]),
                "protein_id": row["protein_id"],
                "start": int(row["start"]),
                "end": int(row["end"])
            }

    return features

####################################################################

#   Annotation loops

####################################################################

def add_signal_peptides_to_genbank(gbk_file, output_file, signal_peptides, lipid_sites, tat_motifs, args):
    records = list(SeqIO.parse(gbk_file, "genbank"))
    count_added = 0

    for record in records:
        for i, feature in enumerate(record.features):
            if feature.type != "CDS":
                continue

            if args.proteinid:
                tags_to_check = feature.qualifiers.get("protein_id", [])
            else:
                tags_to_check = feature.qualifiers.get("locus_tag", []) + feature.qualifiers.get("old_locus_tag", [])
                
            matching_tag = None
            for tag in tags_to_check:
                if tag in signal_peptides or tag in lipid_sites or tag in tat_motifs:
                    matching_tag = tag
                    break
            if not matching_tag:
                continue

            strand = feature.location.strand
            cds_start = int(feature.location.start)
            cds_end = int(feature.location.end)

            # -----------------------
            # Add signal peptide (if any) - calculate amino acid starts and stops
            # -----------------------
            if matching_tag in signal_peptides:
                n_start_aa, h_end_aa = signal_peptides[matching_tag]

                if strand == 1:
                    sig_start = cds_start + (n_start_aa - 1) * 3
                    sig_end = cds_start + (h_end_aa * 3)
                elif strand == -1:
                    sig_end = cds_end - (n_start_aa - 1) * 3
                    sig_start = cds_end - (h_end_aa * 3)
                else:
                    continue

                sig_start = max(sig_start, 0)
                sig_end = min(sig_end, len(record.seq))

                location = FeatureLocation(sig_start, sig_end, strand=strand)
                new_feature = SeqFeature(
                    location=location,
                    type="misc_feature",
                    qualifiers={
                        "note": ["signal peptide"],
                        "locus_tag": [matching_tag]
                    }
                )
                record.features.insert(i + 1, new_feature)
                count_added += 1

            # -------------------------------
            # Add lipid-modified cysteine (if any)
            # -------------------------------
            if matching_tag in lipid_sites:
                aa_pos = lipid_sites[matching_tag]

                if strand == 1:
                    lip_start = cds_start + (aa_pos - 1) * 3
                    lip_end = lip_start + 3
                elif strand == -1:
                    lip_end = cds_end - (aa_pos - 1) * 3
                    lip_start = lip_end - 3
                else:
                    continue

                lip_start = max(lip_start, 0)
                lip_end = min(lip_end, len(record.seq))

                location = FeatureLocation(lip_start, lip_end, strand=strand)
                lipid_feature = SeqFeature(
                    location=location,
                    type="misc_feature",
                    qualifiers={
                        "note": ["lipid-modified cysteine"],
                        "locus_tag": [matching_tag]
                    }
                )
                record.features.insert(i + 1, lipid_feature)
                count_added += 1

            # -------------------------------
            # Add twin-arginine motif (if any)
            # -------------------------------
            if matching_tag in tat_motifs:
                aa_start, aa_end = tat_motifs[matching_tag]

                if strand == 1:
                    tat_start = cds_start + (aa_start - 1) * 3
                    tat_end = cds_start + (aa_end * 3)
                elif strand == -1:
                    tat_end = cds_end - (aa_start - 1) * 3
                    tat_start = cds_end - (aa_end * 3)
                else:
                    continue

                tat_start = max(tat_start, 0)
                tat_end = min(tat_end, len(record.seq))

                location = FeatureLocation(tat_start, tat_end, strand=strand)
                tat_feature = SeqFeature(
                    location=location,
                    type="misc_feature",
                    qualifiers={
                        "note": ["twin-arginine motif"],
                        "locus_tag": [matching_tag]
                    }
                )
                record.features.insert(i + 1, tat_feature)
                count_added += 1

    with open(output_file, "w") as out_handle:
        SeqIO.write(records, out_handle, "genbank")

    print(f"✅ Added annotations to {count_added} CDS features.")
    
###################################################################   

def add_beta_barrels_to_genbank(gbk_file, output_file, beta_features, args):
    """
    Add a single misc_feature per CDS with >=8 Beta sheets.
    The feature spans from first to last Beta sheet residue and notes the strand count.
    """
    records = list(SeqIO.parse(gbk_file, "genbank"))
    count_annotated = 0

    for record in records:
        for feature in record.features:
            if feature.type != "CDS":
                continue

            if args.proteinid:
                feature_keys = feature.qualifiers.get("protein_id", [])
            else:
                feature_keys = feature.qualifiers.get("locus_tag", []) + feature.qualifiers.get("old_locus_tag", [])

            matching_tag = next((tag for tag in feature_keys if tag in beta_features), None)
            
            if not matching_tag:
                continue
            
            region_list = beta_features[matching_tag]
            beta_sheets = [r for r in region_list if r["type"].strip().lower() == "beta sheet"]
            
            if len(beta_sheets) < 4:
                continue  # skip if fewer than 8 beta sheets

            # print(f"✅ Annotating {matching_tag} with {len(beta_sheets)} beta strands")
            # Get residue start/end from first and last beta sheet
            aa_start = beta_sheets[0]["start"]
            aa_end = beta_sheets[-1]["end"]

            # Convert to nucleotide coordinates
            cds_start = int(feature.location.start)
            strand = feature.location.strand
            protein_length = len(feature.qualifiers.get("translation", [""])[0])

            if strand == 1:
                nuc_start = cds_start + (aa_start - 1) * 3
                nuc_end = cds_start + (aa_end) * 3
            else:
                # For reverse strand, reverse AA offset from CDS start
                nuc_end = cds_start + (protein_length - (aa_start - 1)) * 3
                nuc_start = cds_start + (protein_length - aa_end) * 3

            location = FeatureLocation(nuc_start, nuc_end, strand=strand)
            note_text = f"{len(beta_sheets)} strand beta barrel"
            new_feature = SeqFeature(
                location=location,
                type="misc_feature",
                qualifiers={"note": [note_text]}
            )
            cds_index = record.features.index(feature)
            record.features.insert(cds_index + 1, new_feature)
            count_annotated += 1

    SeqIO.write(records, output_file, "genbank")
    print(f"✅ Added beta barrel annotations to {count_annotated} CDS features.")

##################################################################   

def add_multihemes_to_genbank(gbk_file, output_file, multihemes):
    """
    Adds a misc_feature to the GenBank file for each multiheme cytochrome entry.

    Args:
        gbk_file: input GenBank file
        output_file: path to write updated GenBank file
        multihemes: dict from parse_multiheme_tsv(), keyed by locus_tag
    """
    records = list(SeqIO.parse(gbk_file, "genbank"))
    count = 0

    for record in records:
        for feature in record.features:
            if feature.type != "CDS":
                continue

            tags = feature.qualifiers.get("locus_tag", []) + feature.qualifiers.get("old_locus_tag", [])
            matching_tag = next((tag for tag in tags if tag in multihemes), None)
            if not matching_tag:
                continue

            data = multihemes[matching_tag]
            strand = feature.location.strand
            location = FeatureLocation(data["start"], data["end"], strand=strand)

            note = f"{data['motifs_total']}-heme cytochrome"
            qualifiers = {
                "note": [note],
                "locus_tag": [matching_tag]
            }

            new_feature = SeqFeature(location=location, type="misc_feature", qualifiers=qualifiers)

            # Insert immediately after the CDS feature
            cds_index = record.features.index(feature)
            seq_len = len(record.seq)
            if not (0 <= data["start"] < data["end"] <= seq_len):
                print(f"⚠️  Skipping {matching_tag}: region {data['start']}..{data['end']} out of bounds for record {record.id} (length {seq_len})")
                continue  
            record.features.insert(cds_index + 1, new_feature)

            # print(f"✅ Annotated {matching_tag} with multiheme feature")
            count += 1

    SeqIO.write(records, output_file, "genbank")
    print(f"✅ Added {count} multiheme cytochrome annotations to GenBank file: {output_file}")

##################################################################   

#   Main loop

####################################################################
# 
# If two folders are provided with features and gbk files, will do batch mode, otherwise will attempt single mode

def main():
    args = parse_arguments()

    ########## Batch mode: if both folders provided
    
    if args.genbank_directory and args.feature_directory:
        gbk_dir = Path(args.genbank_directory)
        feature_dir = Path(args.feature_directory)
    
        print(f"📁 Batch mode: scanning {feature_dir} for feature files...")
    
        feature_files = list(feature_dir.glob("*.gff3" if args.signal or args.beta else "*.tsv"))
    
        if not feature_files:
            sys.exit("❌ No feature files found in feature directory.")
    
        for feature_file in feature_files:
            print(f"📄 Checking feature file: {feature_file.name}")
    
            # Use regex to extract the accession prefix
            match = re.match(r"^(GCF_\d+\.\d+)", feature_file.name)
            if not match:
                print(f"⚠️ Cannot determine accession from: {feature_file.name}")
                continue
    
            prefix = match.group(1)
            print(f"Using prefix: {prefix}")
    
          # Match either exact or suffix-extended GenBank file, like GCF_99999999.gbk and GCF_99999999_signal.gbk
            matching_gbk = next(
                iter(sorted(gbk_dir.glob(f"{prefix}.gb*")) + sorted(gbk_dir.glob(f"{prefix}_*.gb*"))),
                None
            )
    
            # Build output filename
            suffix = "_signal.gbk" if args.signal else "_beta.gbk" if args.beta else "_multiheme.gbk"
            output_path = gbk_dir / f"{prefix}{suffix}"
    
            # Skip if exists and not forcing overwrite
            if output_path.exists() and not args.force:
                print(f"⏩ Skipping {prefix}: output file {output_path.name} already exists (use --force to overwrite)")
                continue
    
            # Dispatch annotation based on mode
            if args.signal:
                annotations = parse_signalp_gff(feature_file)
                signal_peptides = extract_signal_peptides(annotations)
                lipid_sites = extract_lipid_cysteines(annotations)
                tat_motifs = extract_tat_motifs(annotations)
                add_signal_peptides_to_genbank(matching_gbk, output_path, signal_peptides, lipid_sites, tat_motifs, args)
    
            elif args.beta:
                beta_barrels = parse_beta_gff(feature_file)
                add_beta_barrels_to_genbank(matching_gbk, output_path, beta_barrels, args)
    
            elif args.multiheme:
                multihemes = parse_multiheme_tsv(feature_file)
                add_multihemes_to_genbank(matching_gbk, output_path, multihemes)
    
        print("✅ Batch annotation complete.")
        return  # prevent single-file mode from running

    ########## Single-file fallback mode
    
    check_input_files(args.gff_file, args.genbank_file, args)

    print(f"✅ Input file: {args.gff_file}")
    print(f"✅ GenBank input file: {args.genbank_file}")
    print(f"➡️ Output will be written to: {args.output}")

    # Dispatch based on analysis type
    if args.signal:
        print("🔍 Parsing SignalP GFF...")
        annotations = parse_signalp_gff(args.gff_file)

        signal_peptides = extract_signal_peptides(annotations)
        lipid_sites = extract_lipid_cysteines(annotations)
        tat_motifs = extract_tat_motifs(annotations)

        add_signal_peptides_to_genbank(
            args.genbank_file,
            args.output,
            signal_peptides,
            lipid_sites,
            tat_motifs
        )

    elif args.beta:
        print("🔍 Parsing beta-barrel GFF...")
        beta_barrels = parse_beta_gff(args.gff_file)

        add_beta_barrels_to_genbank(
            args.genbank_file,
            args.output,
            beta_barrels
        )

    elif args.multiheme:
        print("🔍 Parsing multiheme cytochrome file...")
        multihemes = parse_multiheme_tsv(args.gff_file)

        add_multihemes_to_genbank(
            args.genbank_file,
            args.output,
            multihemes
        )

    print(f"✅ All annotations added to GenBank file: {args.output}")
    
if __name__ == "__main__":
    main()