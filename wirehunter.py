#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wirehunter.py

Identify gene clusters in a GenBank genome based on homology to user-supplied protein sequences.

Input:
- GenBank genome file (.gbk)
- Protein query file in TSV format with columns: ID, name, sequence

Output:
- GenBank file(s) of extracted regions
- Optional: intermediate BLAST database/results

By default, temporary BLAST files are deleted after use unless --keep-temp is specified.

Example usage:
    python wirehunter.py genome.gbk queries.tsv
    
Example to alter the number of hits that makes a postive, and window to find the genes 

    wirehunter.py genome.gbk queries.tsv --min_hits=3 --window=15000
"""


import os
import sys
import argparse
import tempfile
import shutil
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import csv



sys.argv = [
    "wirehunter.py",
    "/Users/daniel/Desktop/ReturnoftheChrome/wirehunter/GCF_000016745.1_multiheme.gbk",
    "/Users/daniel/Desktop/ReturnoftheChrome/wirehunter/omce.tsv",
    "-o", "/Users/daniel/Desktop/ReturnoftheChrome/wirehunter/GCF_000016745.1_results",
    "--window=4",
    "--min_hits=3",
]



def parse_args():
    parser = argparse.ArgumentParser(description="Search for gene clusters based on homologous proteins")
    parser.add_argument("genbank_file", help="Path to GenBank genome file")
    parser.add_argument("query_file", help="TSV file with protein sequences to search")
    parser.add_argument("--output", "-o", default="cluster_results", help="Output directory for results")
    parser.add_argument("--window", type=int, default=2, help="Number of genes to include upstream/downstream of cluster")
    parser.add_argument("--min_hits", type=int, default=4, help="Minimum distinct query IDs needed to define a cluster")
    parser.add_argument("--max-distance", type=int, default=10000, help="Max genomic distance (bp) allowed between hits in a cluster")
    parser.add_argument("--keep-temp", action="store_true", help="Keep intermediate BLAST database and results")
    return parser.parse_args()

def read_query_tsv(tsv_file):
    queries = {}
    with open(tsv_file) as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        if header != ['ID', 'name', 'sequence']:
            print(f"⚠️ Warning: Header is not exactly 'ID\\tname\\tsequence'. Using positions instead.")

        for line_num, row in enumerate(reader, start=2):
            if len(row) < 3:
                print(f"⚠️ Skipping malformed line {line_num}: {row}")
                continue
            group_id = row[0].strip()
            name = row[1].strip()
            seq = row[2].replace(" ", "").strip()
            if not group_id or not name or not seq:
                print(f"⚠️ Incomplete entry on line {line_num}")
                continue
            queries.setdefault(group_id, []).append((name, seq))
    return queries

def extract_proteins_from_genbank(gbk_file):
    proteins = []
    records = list(SeqIO.parse(gbk_file, "genbank"))
    for record in records:
        for feature in record.features:
            if feature.type == "CDS" and "translation" in feature.qualifiers:
                locus = feature.qualifiers.get("locus_tag", [None])[0]
                if not locus:
                    continue
                translation = feature.qualifiers["translation"][0]
                # Explicitly format FASTA header
                record = SeqRecord(
                    Seq(translation),
                    id=locus,
                    description=locus  # required to prevent gnl|BL_ORD_ID| overwrite
                )
                proteins.append(record)
    return proteins

def write_fasta(seqs, path):
    with open(path, "w") as out:
        for record in seqs:
            header = record.id.strip()  # should be locus_tag like GSU_1234
            out.write(f">{header}\n{str(record.seq)}\n")
            
def run_blastp(query_fasta, db_path, temp_dir):
    out_xml = os.path.join(temp_dir, f"{Path(query_fasta).stem}_blast.xml")
    blast_cline = NcbiblastpCommandline(query=query_fasta, db=db_path, evalue=1e-5,
                                        outfmt=5, out=out_xml)
    stdout, stderr = blast_cline()
    return out_xml

def parse_blast_results(xml_path, group_id=None, group_name=None):
    """
    Parse BLAST XML results and return a dict of hits by query.
    Also prints the best hit for the group, if found.

    Returns:
        dict: query_id -> list of subject locus_tags (cleaned)
    """
    hits_by_query = {}
    best_hit = None
    best_identity = 0.0
    best_subject = ""

    with open(xml_path) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
    
        for record in blast_records:
            qid = record.query.split()[0]
            best_hsp = None
            best_score = -1
            best_subject = None
            best_identity = 0.0
    
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    if hsp.score > best_score:
                        best_score = hsp.score
                        best_subject = alignment.hit_id
                        best_identity = (hsp.identities / hsp.align_length) * 100
    
            if best_subject:
                hits_by_query[qid] = [best_subject]
                print(f" Best hit to {group_name} (ID {group_id}): {best_subject} ({best_identity:.1f}% identity)")
            else:
                print(f" ! No hits found for {group_name} (ID {group_id})")

    # Extract locus-like part of the best subject
    # Keep full subject ID
        if best_subject:
            print(f" Best hit to {group_name} (ID {group_id}): {best_subject} ({best_identity:.1f}% identity)")
        else:
            print(f" ! No hits found for {group_name} (ID {group_id})")
        
        # Return IDs as-is (they are already clean locus_tags)
        return hits_by_query

def extract_region_from_gbk(gbk_path, loci_to_extract, window):
    out_records = []
    for record in SeqIO.parse(gbk_path, "genbank"):
        loci_data = []
        for i, feature in enumerate(record.features):
            if feature.type != "CDS":
                continue
            locus = feature.qualifiers.get("locus_tag", [None])[0]
            if locus:
                start = int(feature.location.start)
                end = int(feature.location.end)
                loci_data.append((locus, start, end, i))

        # Map locus to index in loci_data
        loci_to_indices = {locus: i for i, (locus, *_rest) in enumerate(loci_data)}

        included_indices = [loci_to_indices[locus] for locus in loci_to_extract if locus in loci_to_indices]

        included_indices = sorted(set(included_indices))
        if not included_indices:
            continue

        # Get span of region from earliest to latest hit gene
        span_start = loci_data[included_indices[0]][1]
        span_end = loci_data[included_indices[-1]][2]
        # print(f"📏 Cluster genomic span: {span_start}–{span_end}")

        # Create new sliced record
        new_record = record[span_start:span_end]
        new_record.id = record.id
        new_record.name = record.name
        new_record.description = record.description
        new_record.annotations = record.annotations.copy()
        if "molecule_type" not in new_record.annotations:
            new_record.annotations["molecule_type"] = "DNA"

        # Include all features in this span, with coordinates shifted
        new_features = []
        for feature in record.features:
            if feature.location.start >= span_start and feature.location.end <= span_end:
                shifted = feature._shift(-span_start)
                new_features.append(shifted)
        new_record.features = new_features
        new_record.annotations["topology"] = "linear"
        out_records.append(new_record)

    return out_records

def find_flexible_cluster(hits_by_group, min_hits, gbk_path, max_distance=10000, window=2):
    from itertools import combinations, product
    from collections import defaultdict
    from Bio import SeqIO

    group_to_loci = defaultdict(set)
    for group_id, loci in hits_by_group.items():
        for locus in loci:
            group_to_loci[group_id].add(locus)

    group_ids = list(group_to_loci.keys())
    if len(group_ids) < min_hits:
        return set()

    # Parse genome CDS features
    loci_data = []
    for record in SeqIO.parse(gbk_path, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            locus = feature.qualifiers.get("locus_tag", [None])[0]
            if not locus:
                continue
            start = int(feature.location.start)
            end = int(feature.location.end)
            loci_data.append((locus, start, end, feature))

    locus_to_info = {locus: (start, end) for locus, start, end, _ in loci_data}

    for combo_size in range(min_hits, len(group_ids) + 1):
        for group_subset in combinations(group_ids, combo_size):
            hit_sets = [group_to_loci[gid] for gid in group_subset]
            for candidate in product(*hit_sets):
                if len(set(candidate)) < min_hits:
                    continue

                coords = [locus_to_info[locus] for locus in candidate if locus in locus_to_info]
                if not coords:
                    continue

                span_start = min(pos[0] for pos in coords)
                span_end = max(pos[1] for pos in coords)
                span = span_end - span_start
                if span > max_distance:
                    continue

                # Extract the exact region by genomic coordinate range
                in_region_indices = [i for i, (_, start, end, _) in enumerate(loci_data)
                                     if start >= span_start and end <= span_end]

                if not in_region_indices:
                    continue

                leftmost_idx = max(0, min(in_region_indices) - window)
                rightmost_idx = min(len(loci_data), max(in_region_indices) + window + 1)

                region_loci = [loci_data[i][0] for i in range(leftmost_idx, rightmost_idx)]

                print(f"🧪 Upstream CDS: {min(in_region_indices) - leftmost_idx} / "
                      f"Downstream CDS: {rightmost_idx - max(in_region_indices) - 1}")
                print(f"Extracting {len(region_loci)} total genes...")
                print(f"<--> Cluster genomic span: {loci_data[leftmost_idx][1]}–{loci_data[rightmost_idx - 1][2]}")

                return set(region_loci)

    return set()

def main():
    args = parse_args()
    os.makedirs(args.output, exist_ok=True)
    temp_dir = tempfile.mkdtemp(dir=args.output)

    print("Reading query proteins...")
    query_groups = read_query_tsv(args.query_file)

    print("Extracting proteins from GenBank file...")
    proteins = extract_proteins_from_genbank(args.genbank_file)
    prot_fasta = os.path.join(temp_dir, "genome_proteins.faa")
    write_fasta(proteins, prot_fasta)

    print("Running BLASTP for each query group...")# Step: Make BLAST DB once
    db_path = os.path.join(temp_dir, "blastdb")
    print(f"Creating BLAST database for genome proteins...")
    subprocess.run([
    "makeblastdb", "-in", prot_fasta, "-dbtype", "prot",
    "-parse_seqids",  # This is key to parsing the sequence IDs
    "-out", db_path
    ], check=True)
    
    all_hits = {}
    for group_id, entries in query_groups.items():
        group_fasta = os.path.join(temp_dir, f"query_{group_id}.faa")
        group_seqs = [SeqRecord(Seq(seq), id=f"{group_id}_{i}", description=name)
    for i, (name, seq) in enumerate(entries)]
        write_fasta(group_seqs, group_fasta)
        result_xml = run_blastp(group_fasta, db_path, temp_dir)
        # Determine a representative name for this group
        group_name = entries[0][0] if entries else f"group_{group_id}"
        hits = parse_blast_results(result_xml, group_id, group_name)
        all_hits[group_id] = [hit for hits in hits.values() for hit in hits]

    print(" Identifying clusters...")
    cluster_loci = find_flexible_cluster(
    hits_by_group=all_hits,
    min_hits=args.min_hits,
    gbk_path=args.genbank_file,
    max_distance=10000,
    window=args.window
    )

    # print(f"📦 Extracting {len(cluster_loci)} context genes...")
    extracted = extract_region_from_gbk(args.genbank_file, cluster_loci, args.window)

    output_gbk = os.path.join(args.output, "cluster_regions.gbk")
    with open(output_gbk, "w") as out:
        SeqIO.write(extracted, out, "genbank")

    if args.keep_temp:
        print(f" Temporary files kept in: {temp_dir}")
    else:
        shutil.rmtree(temp_dir)
        print(" Temporary files cleaned up.")

    print(f"✅ Done. Extracted regions written to {output_gbk}")

if __name__ == "__main__":
    main()
