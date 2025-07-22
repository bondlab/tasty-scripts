#!/usr/bin/env python3

import argparse
import os
import re
import sys
import csv
# import glob
from pathlib import Path
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

#
# Featurebot - read GFF/TSV/TXT files resulting from various prediction tools
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
#   Next job for featurebot is to process TSV best hit data and label genes
#
#
# Simulate command-line arguments when running in an environment like Spyder.
# comment out for command line operation


# sys.argv = [
#     "featurebot.py",
#     "-gbkdir", "/Users/daniel/Projects/wirehunter/all_3354_delta_wirehunter_results/1_cluster_annotated_bgk/",
#     "-featuredir", "/Users/daniel/Projects/wirehunter/all_3354_delta_wirehunter_results/2_heme_tsvs/",
#     "-multiheme",
# #    "-proteinid",
# #    "--force"
# ]

# # # /Users/daniel/Projects/finals/TypeI_results_3hits/
# #
####################################################################


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Featurebot: map four flavors of results from analyses on to genome annotations.\n\n"
                    "Matches outputs from SignalP, beta barrel TMHMM, multiheme analyses, or cluster protein finders\n"
                    "with their respective GenBank files, when both have filenames beginning with the same accession.\n\n"
                    "Invoke batch mode using -gbkdir abnd -featuredir, single analyses with -gbk and -gff.\n"
                    "Invoke different analyses with -signal, -beta, -multiheme, or -cluster.\n"
                    "Default is to look for a locus tag for matching, switch to proteinID using -proteinid.\n"
                    "Features are mapped as misc_feature with coordinates when available.\n"
    )

    parser.add_argument("-gff", "--gff_file", required=False,
                        help="Path to feature file (GFF3/TSV/cluster TSV).")
    parser.add_argument("-gbk", "--genbank_file", required=False,
                        help="Path to GenBank (.gb/.gbk/.gbff) file.")
    parser.add_argument("-gbkdir", "--genbank_directory", required=False,
                        help="Directory with GenBank files for batch mode.")
    parser.add_argument("-featuredir", "--feature_directory", required=False,
                        help="Directory with feature files for batch mode.")
    parser.add_argument("--force", action="store_true",
                        help="Force overwrite of existing output files.")
    parser.add_argument("-proteinid", action="store_true",
                        help="Use protein_id instead of locus_tag for matching.")


    parser.add_argument("--suffix", type=str, default=None,
                    help="Suffix to add to output filenames and output directory (e.g., 'signal', 'beta', 'cluster'). Overrides default suffix.")
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-signal", action="store_true", help="Parse SignalP predictions.")
    group.add_argument("-beta", action="store_true", help="Parse beta-barrel predictions.")
    group.add_argument("-multiheme", action="store_true", help="Parse multiheme predictions.")
    group.add_argument("-cluster", action="store_true", help="Parse cluster protein mapping data.")

    
    args = parser.parse_args()

    # Auto-generate output file name in single mode
    # base, _ = os.path.splitext(os.path.basename(args.genbank_file))
    # suffix = "_signal.gbk" if args.signal else "_beta.gbk" if args.beta else "_multiheme.gbk" if args.multiheme else "_cluster.gbk"
    # output_file = base + suffix

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
        if not os.path.isfile(gbk_path) or not gbk_path.lower().endswith((".gff", ".gff3")):
            sys.exit(f"❌ DEBUG Error: expected a .gff or .gff3 file for signal or beta analysis, got: {gff_path}")

    if not os.path.isfile(gbk_path) or not gbk_path.lower().endswith((".gb", ".gbk", "gbff")):
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

# cluster output parsing ###################################################################

def parse_cluster_tsv(cluster_file):
    cluster_annotations = {}
    with open(cluster_file, newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
#         header = next(reader, None)
        for row in reader:
            if len(row) < 3:
                continue
            locus_tag, protein_name, organism = row[0].strip(), row[1].strip(), row[2].strip()

            # Skip entries where protein name is "No hit"
            if protein_name.lower() == "no hit":
                continue

            # Skip incomplete rows
            if not locus_tag or not protein_name or not organism:
                continue

            # Save protein_name and organism as a combined string
            cluster_annotations[locus_tag] = f"{protein_name} {organism}"

    return cluster_annotations

####################################################################

#   Annotation loops - add features to GenBank, as misc_feature and/or /product (clusters)

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
                continue  # skip if fewer than n beta sheets

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

def add_clusters_to_genbank(gbk_file, output_file, cluster_data):
    records = list(SeqIO.parse(gbk_file, "genbank"))
    count = 0
    for record in records:
        for feature in record.features:
            if feature.type != "CDS":
                continue
            locus_tags = feature.qualifiers.get("locus_tag", []) + feature.qualifiers.get("old_locus_tag", [])
            matching_tag = next((tag for tag in locus_tags if tag in cluster_data), None)
            if not matching_tag:
                continue
            combined_note = cluster_data[matching_tag]

            #
            #
            #
            # Replace /product with the new label (can comment this out if not needed)
            #
            feature.qualifiers["product"] = [combined_note]
            # 


            # Add as a misc_feature as well
            new_feature = SeqFeature(
                location=feature.location,
                type="misc_feature",
                qualifiers={
                    "note": [combined_note],
                    "locus_tag": [matching_tag]
                }
            )
            cds_index = record.features.index(feature)
            record.features.insert(cds_index + 1, new_feature)
            count += 1
    with open(output_file, "w") as out_handle:
        SeqIO.write(records, out_handle, "genbank")
        
    print(f"Added {count} cluster annotations (and updated /product) to GenBank file: {output_file}")
  
##################################################################

#   Main loop

####################################################################

# If two folders are provided with features and gbk files, will do batch mode, otherwise will attempt single mode


def main():
    args = parse_arguments()

    # Determine suffix from argument or mode
    if args.suffix:
        suffix = f"_{args.suffix}"
    elif args.signal:
        suffix = "_signal"
    elif args.beta:
        suffix = "_beta"
    elif args.multiheme:
        suffix = "_multiheme"
    elif args.cluster:
        suffix = "_cluster"
    else:
        suffix = "_annotated"

    ########## Batch mode: if both folders provided
    if args.genbank_directory and args.feature_directory:
        gbk_dir = Path(args.genbank_directory)
        feature_dir = Path(args.feature_directory)
    
        print(f"📁 Batch mode: scanning {feature_dir} for feature files...")
    
        # Match GFF or TSV depending on mode
        feature_files = list(
            feature_dir.glob("*.gff3" if args.signal or args.beta else "*.tsv")
        )
        if not feature_files:
            sys.exit("❌ No feature files found in feature directory.")
    
        # Create output directory next to GenBank directory
        outdir = gbk_dir.parent / f"{gbk_dir.name}{suffix}"
        outdir.mkdir(parents=True, exist_ok=True)
    
        for feature_file in feature_files:
            print(f"\n Checking feature file: {feature_file.name}")
    
            # Extract accession (e.g., GCA_000016745.1)
            match = re.match(r"^(GCF|GCA)_\d+\.\d+", feature_file.name)
            if not match:
                print(f"⚠️ Cannot extract accession from: {feature_file.name}")
                continue
            accession = match.group(0)
            print(f" Using accession: {accession}")
    
            # Find GenBank file starting with accession
            gbk_matches = sorted(gbk_dir.glob(f"{accession}*.gb*"))
            if not gbk_matches:
                print(f"⚠️ No GenBank file found starting with {accession}")
                print(" GenBank files available:")
                for f in gbk_dir.glob("*.gb*"):
                    print(f"  - {f.name}")
                continue
            matching_gbk = gbk_matches[0]
            print(f"✅ Matched GenBank file: {matching_gbk.name}")
    
            # Build output path: accession + suffix
            output_path = outdir / f"{accession}{suffix}.gbk"
    
            if output_path.exists() and not args.force:
                print(f"⏩ Skipping {accession}: output file {output_path.name} already exists (use --force to overwrite)")
                continue
    
            # Dispatch by mode
            if args.signal:
                annotations = parse_signalp_gff(feature_file)
                signal_peptides = extract_signal_peptides(annotations)
                lipid_sites = extract_lipid_cysteines(annotations)
                tat_motifs = extract_tat_motifs(annotations)
                add_signal_peptides_to_genbank(
                    matching_gbk, output_path, signal_peptides, lipid_sites, tat_motifs, args
                )
            elif args.beta:
                beta_barrels = parse_beta_gff(feature_file)
                add_beta_barrels_to_genbank(matching_gbk, output_path, beta_barrels, args)
            elif args.multiheme:
                multihemes = parse_multiheme_tsv(feature_file)
                add_multihemes_to_genbank(matching_gbk, output_path, multihemes)
            elif args.cluster:
                cluster_data = parse_cluster_tsv(feature_file)
                add_clusters_to_genbank(matching_gbk, output_path, cluster_data)
    
        print(f"\n Batch annotation complete. Output files written to: {outdir}")
        return
    ########## Single-file mode fallback

    # check_input_files(args.gff_file, args.genbank_file, args)

    gbk_path = Path(args.genbank_file)
    base = gbk_path.stem
    parent_dir = gbk_path.parent

    outdir = parent_dir / f"{parent_dir.name}{suffix}"
    outdir.mkdir(parents=True, exist_ok=True)

    output_file = outdir / f"{base}{suffix}.gbk"

    print(f"Input file: {args.gff_file}")
    print(f"GenBank input file: {args.genbank_file}")
    print(f"➡️ Output will be written to: {output_file}")

    if args.signal:
        print("🔍 Parsing SignalP GFF...")
        annotations = parse_signalp_gff(args.gff_file)
        signal_peptides = extract_signal_peptides(annotations)
        lipid_sites = extract_lipid_cysteines(annotations)
        tat_motifs = extract_tat_motifs(annotations)
        add_signal_peptides_to_genbank(
            args.genbank_file, output_file, signal_peptides, lipid_sites, tat_motifs, args
        )
    elif args.beta:
        print("🔍 Parsing beta-barrel GFF...")
        beta_barrels = parse_beta_gff(args.gff_file)
        add_beta_barrels_to_genbank(
            args.genbank_file, output_file, beta_barrels, args
        )
    elif args.multiheme:
        print("🔍 Parsing multiheme cytochrome file...")
        multihemes = parse_multiheme_tsv(args.gff_file)
        add_multihemes_to_genbank(
            args.genbank_file, output_file, multihemes
        )
    elif args.cluster:
        print("🔍 Parsing cluster TSV...")
        cluster_data = parse_cluster_tsv(args.gff_file)
        add_clusters_to_genbank(
            args.genbank_file, output_file, cluster_data
        )

    print(f"All annotations added to GenBank file: {output_file}")

##################################################################
if __name__ == "__main__":
    main()
    