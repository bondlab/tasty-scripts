#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wirehunter_0.5.py

The improved version that runs in batch mode too!


Identify gene clusters in GenBank genomes based on homology to user-supplied protein sequences.
Supports:
- Single GenBank file input
- Batch folder mode (--batch)
- Cluster summary table logging accession, LOCUS, ORGANISM, DEFINITION, and cluster presence (0/1)

Output:
- GenBank file(s) of extracted regions
- Table of best hits to each protein in the cluster
    These two files can be combined by the Featurebot program to annotate the hits
    
- Logfile of all hits to the query database
- Cluster summary TSV across all processed genomes

- Optional: intermediate BLAST database/results
By default, temporary BLAST files are deleted after use unless --keep-temp is specified.

Example usage: basic level
    python wirehunter.py genome.gbk queries.tsv
    
Example with adjusted window and minimum hits:
    python wirehunter.py genome.gbk queries.tsv --min_hits=3 --window=15000
    
Example pushing all the buttons: require hits to 4 proteins in the database, with 80% coverage, all hits must be within 15000 bp contig
    wirehunter.py database.tsv \
    -o results_4_hits_70 \
    --batch genomes \
    --window=4 --max-distance=15000 --min_hits=4 --min-coverage=80


It will only return putative clusters that are on the same contig. 

Example output table (cluster_best_hits.tsv):
locus_tag    best_query_name    organism               protein_id       identity
GSU_1234     OmcE-like          G. metallireducens     WP_123456        72.4


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
from collections import defaultdict
from itertools import combinations, product
import csv

# -----------------------------------------
# Comment out the block below for normal usage
# use for quick debugging in IDEs like Spyder

#    Change Protein ID of found OmcE to Locus tag - eaisier to see if 
#       you're in the same region
#
#   Comment out ! No hits found 
#   and Hit for updates, save only Group 1 --> in output



# sys.argv = [
#     "wirehunter_0.5.py",
#     "/Users/daniel/Projects/wirehunter/OmcE_database_Type_II.tsv",
#     "-o", "/Users/daniel/Projects/Downloaded_genomes/Deltas_RefSeq_848/TypeII_4_update",
#     "--batch",
#     "/Users/daniel/Projects/Downloaded_genomes/Deltas_RefSeq_848/Delta_gbff/",
#     "--window=4",
#     "--max-distance=15000",
#     "--min_hits=4",
#     "--min-coverage=70"
# ]

# -----------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description="Search for gene clusters based on homologous proteins")
    parser.add_argument("query_file", help="5-column TSV file with protein sequences to search")
    parser.add_argument("--genbank_file", help="Path to a GenBank genome file")
    parser.add_argument("--batch", help="Folder containing GenBank files to process in batch mode")
    parser.add_argument("--output", "-o", default="cluster_results", help="Output directory for results")
    parser.add_argument("--window", "-w", type=int, default=2, help="Number of genes to include upstream/downstream of cluster")
    parser.add_argument("--min_hits", "-mh", type=int, default=4, help="Minimum distinct query IDs needed to define a cluster")
    parser.add_argument("--max-distance", type=int, default=10000, help="Max genomic distance (bp) allowed between hits in a cluster")
    parser.add_argument("--keep-temp", action="store_true", help="Keep intermediate BLAST database and results")
    parser.add_argument("--min-coverage", type=float, default=70.0, help="Minimum coverage (%%) for query and subject proteins")
    parser.add_argument("--min-identity", type=float, default=35.0, help="Minimum identity (%%) for BLAST hits")
    return parser.parse_args()


def read_query_tsv(tsv_file):
    """
    Read a TSV file with 5 columns:
    ID, name, organism, protein_id, sequence
    Returns:
        dict: group_id -> list of (name, seq, organism, protein_id)
    """
    queries = {}
    with open(tsv_file) as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        expected_header = ['ID', 'name', 'organism', 'protein_id', 'sequence']
        if header != expected_header:
            print(f"⚠️ Warning: Expected header {expected_header}, got {header}. Will use positions instead.")

        for line_num, row in enumerate(reader, start=2):
            if len(row) < 5:
                print(f"⚠️ Skipping malformed line {line_num}: {row}")
                continue
            group_id, name, organism, protein_id, seq = [x.strip() for x in row[:5]]
            seq = seq.replace(" ", "")
            if not group_id or not name or not seq:
                print(f"⚠️ Incomplete entry on line {line_num}")
                continue
            queries.setdefault(group_id, []).append((name, seq, organism, protein_id))
    return queries

def extract_proteins_from_genbank(gbk_file):
    """
    Extract all CDS proteins from a GenBank file.
    Returns a list of SeqRecord objects.
    """
    proteins = []
    records = list(SeqIO.parse(gbk_file, "genbank"))
    for record in records:
        for feature in record.features:
            if feature.type == "CDS" and "translation" in feature.qualifiers:
                locus = feature.qualifiers.get("locus_tag", [None])[0]
                if not locus or not feature.location:
                    continue
                translation = feature.qualifiers["translation"][0]
                proteins.append(
                    SeqRecord(Seq(translation), id=locus, description=locus)
                )
    return proteins

def write_fasta(seqs, path):
    """
    Write a list of SeqRecords to a FASTA file, ensuring unique IDs.
    Appends _2, _3, etc. to duplicate IDs, prevents script crashing
    from bad GenBank files.
    """
    id_counts = defaultdict(int)
    with open(path, "w") as out:
        for record in seqs:
            base_id = record.id.strip()
            id_counts[base_id] += 1
            if id_counts[base_id] > 1:
                unique_id = f"{base_id}_{id_counts[base_id]}"
            else:
                unique_id = base_id
            out.write(f">{unique_id}\n{str(record.seq)}\n")

def run_blastp(query_fasta, db_path, temp_dir):
    out_xml = Path(temp_dir) / f"{Path(query_fasta).stem}_blast.xml"
    blast_cline = NcbiblastpCommandline(
        query=str(query_fasta), db=str(db_path), evalue=1e-5,
        outfmt=5, out=str(out_xml)
    )
    try:
        stdout, stderr = blast_cline()
    except Exception as e:
        print(f"⚠️ BLAST execution failed for {query_fasta}: {e}")
        sys.exit(1)
    return out_xml

def extract_region_from_gbk(gbk_path, loci_to_extract, window):
    """
    Extract a GenBank region surrounding the given loci, expanded by `window` genes.
    Supports circular genomes (wrap-around handling).
    """
    out_records = []
    for record in SeqIO.parse(gbk_path, "genbank"):
        is_circular = record.annotations.get("topology", "").lower() == "circular"
        seq_len = len(record.seq)

        # Collect CDS loci
        loci_data = []
        for i, feature in enumerate(record.features):
            if feature.type != "CDS":
                continue
            locus = feature.qualifiers.get("locus_tag", [None])[0]
            if not locus or not feature.location:
                continue
            start, end = int(feature.location.start), int(feature.location.end)
            loci_data.append((locus, start, end, i))

        loci_to_indices = {locus: i for i, (locus, *_rest) in enumerate(loci_data)}
        included_indices = sorted({loci_to_indices[locus] for locus in loci_to_extract if locus in loci_to_indices})
        
        if not included_indices:
            continue

        left_idx = min(included_indices) - window
        right_idx = max(included_indices) + window + 1

        if is_circular:
            # Wrap indices around for circular genomes
            wrapped_indices = [(i % len(loci_data)) for i in range(left_idx, right_idx)]
            wrapped_indices = sorted(set(wrapped_indices))
        else:
            # Clamp indices for linear genomes
            left_idx = max(0, left_idx)
            right_idx = min(len(loci_data), right_idx)
            wrapped_indices = list(range(left_idx, right_idx))
        
        # Calculate the start and end positions for the region
        region_starts = [loci_data[i][1] for i in wrapped_indices]
        region_ends = [loci_data[i][2] for i in wrapped_indices]
        
        # Determine genomic span (linearized view)
        span_start = min(region_starts)
        span_end = max(region_ends)

        print(f"  Extracting {len(wrapped_indices)} genes from {span_start}–{span_end}")

        # Determine genomic span
        span_start = min(region_starts)
        span_end = max(region_ends)

        if is_circular and span_end - span_start > seq_len / 2:
            # Region wraps around origin
            span_start, span_end = span_end, span_start  # flipped
            print(" -O- Wrapping around circular genome origin")

            # Concatenate slices: end → start
            new_seq = record.seq[span_start:] + record.seq[:span_end]
        else:
            new_seq = record.seq[span_start:span_end]

        new_record = SeqRecord(new_seq, id=record.id, name=record.name,
                               description=record.description,
                               annotations=record.annotations.copy())
        new_record.annotations["topology"] = "linear"  # extracted region is linear
        new_record.annotations.setdefault("molecule_type", "DNA")

        # Include only features within the extracted region
        new_features = []
        for feature in record.features:
            if not feature.location:
                continue
            f_start, f_end = int(feature.location.start), int(feature.location.end)
            if is_circular and span_end < span_start:
            
            # wrapping region
                if f_start >= span_start or f_end <= span_end:
                    shifted_feature = feature._shift(-span_start % seq_len)
                    new_features.append(shifted_feature)
            else:
                # non-wrapping region
                if f_start >= span_start and f_end <= span_end:
                    shifted_feature = feature._shift(-span_start)
                    new_features.append(shifted_feature)

        new_record.features = new_features
        out_records.append(new_record)
        
        locus_tags_in_region = []
        for i in wrapped_indices:
            if 0 <= i < len(loci_data):
                locus, *_ = loci_data[i]
                locus_tags_in_region.append(locus)

    return out_records, locus_tags_in_region

def find_flexible_cluster(hits_by_group, min_hits, gbk_path, max_distance=10000, window=2):
    """
    Identify clusters containing at least `min_hits` unique query groups.
    Clusters must be on the same contig.
    Adds debug prints to track locus and contig mappings.
    """
    group_to_loci = defaultdict(set)
    locus_to_record = {}  # track which contig each locus is on
    loci_data = []

    # Build locus data: track contig for each locus
    for record in SeqIO.parse(gbk_path, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            locus = feature.qualifiers.get("locus_tag", [None])[0]
            if not locus or not feature.location:
                continue
            start, end = int(feature.location.start), int(feature.location.end)
            loci_data.append((locus, start, end, feature, record.id))
            locus_to_record[locus] = record.id

    locus_to_info = {locus: (start, end, rec_id) for locus, start, end, _, rec_id in loci_data}

    # Debug: show mapping of hits by group
    print("--Hits by group:")
    for group, loci in hits_by_group.items():
        print(f"  Group {group} -> {loci}")

    # # Debug: show locus to contig mapping
    # print("\n🔎 Locus to contig mapping:")
    # for locus in locus_to_info:
    #     start, end, contig = locus_to_info[locus]
    #     print(f"  {locus}: {contig} ({start}-{end})")

    for group_id, loci in hits_by_group.items():
        for locus in loci:
            group_to_loci[group_id].add(locus)

    group_ids = list(group_to_loci.keys())
    if len(group_ids) < min_hits:
        print(f"⚠️ Fewer than {min_hits} genes found: {group_ids}")
        return set()

    for combo_size in range(min_hits, len(group_ids) + 1):
        for group_subset in combinations(group_ids, combo_size):
            hit_sets = [group_to_loci[gid] for gid in group_subset]
            for candidate in product(*hit_sets):
                # Instead of counting unique loci, count unique group IDs represented
                candidate_groups = {
                    group_ids[i] for i, locus in enumerate(candidate) if locus in group_to_loci[group_ids[i]]
                }
                if len(candidate_groups) < min_hits:
                    continue

                # Check that all hits are on the same contig
                contigs = {locus_to_record[locus] for locus in candidate if locus in locus_to_record}
                if len(contigs) > 1:
                    continue  # skip if hits are on multiple contigs

                coords = [locus_to_info[locus][:2] for locus in candidate if locus in locus_to_info]
                if not coords:
                    continue
                span_start = min(pos[0] for pos in coords)
                span_end = max(pos[1] for pos in coords)
                if (span_end - span_start) > max_distance:
                    continue

                # Identify window genes on that contig
                contig_id = next(iter(contigs))
                in_region_indices = [i for i, (locus, start, end, _, rec_id) in enumerate(loci_data)
                                     if rec_id == contig_id and start >= span_start and end <= span_end]
                if not in_region_indices:
                    continue

                left_idx = max(0, min(in_region_indices) - window)
                right_idx = min(len(loci_data), max(in_region_indices) + window + 1)
                region_loci = [loci_data[i][0] for i in range(left_idx, right_idx) if loci_data[i][4] == contig_id]

                print(f"Cluster region includes {len(region_loci)} genes on contig {contig_id}")
                return set(region_loci)

    print("⚠️ No cluster found after checking all combinations.")
    return set()

def log(message, log_handle, also_print=True):
    log_handle.write(message + "\n")
    log_handle.flush()
    if also_print:
        print(message)

def parse_blast_results(xml_path, group_id, group_name, queries, min_coverage=70.0, log_handle=None, also_print=True, min_identity=35.0, locus_tag_map=None):
    """
    Parse BLAST XML results, applying minimum coverage and identity thresholds.
    Tracks the best hit per locus_tag (not just per query!) by bit-score.
    """
    hits_by_query = {}
    best_hits_per_locus = {}
    db_entries = queries.get(group_id, [])
    
    with open(xml_path) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        for record in blast_records:
            qid = record.query.split()[0]

            for alignment in record.alignments:
                subject_id = alignment.hit_id.split()[0]
                locus_tag = locus_tag_map.get(subject_id, subject_id)

                for hsp in alignment.hsps:
                    bits_per_residue = hsp.bits / hsp.align_length
                    query_length = record.query_length
                    subject_length = alignment.length
                    query_coverage = (hsp.align_length / query_length) * 100
                    subject_coverage = (hsp.align_length / subject_length) * 100

                    if query_coverage < min_coverage or subject_coverage < min_coverage:
                        continue

                    identity = (hsp.identities / hsp.align_length) * 100
                    if identity < min_identity:
                        continue

                    # Determine if this is the best hit for this locus_tag
                    if (locus_tag not in best_hits_per_locus) or (hsp.bits > best_hits_per_locus[locus_tag]["bit_score"]):
                        # Identify the associated query protein information
                        if db_entries:
                            try:
                                _, idx_str = qid.rsplit("_", 1)
                                idx = int(idx_str)
                                name, seq, organism, protein_id = db_entries[idx]
                            except (ValueError, IndexError):
                                name = "unknown"
                                organism = "unknown"
                                protein_id = "unknown"
                        else:
                            name, organism, protein_id = "unknown", "unknown", "unknown"

                        best_hits_per_locus[locus_tag] = {
                            "query_id": qid,
                            "query_name": name,
                            "query_organism": organism,
                            "query_protein_id": protein_id,
                            "identity": identity,
                            "bit_score": hsp.bits,
                            "bits_per_residue": bits_per_residue,
                            "query_coverage": query_coverage,
                            "subject_coverage": subject_coverage,
                        }

                        # Logging this as a best hit update
                        log(f" Hit for {locus_tag}: {name} from {organism} "
                            f"({identity:.1f}% identity, query coverage {query_coverage:.1f}, subject coverage {subject_coverage:.1f})",
                            log_handle, also_print)

            # Use record.alignments to ensure no UnboundLocalError
            if record.alignments:
                hits_by_query.setdefault(qid, []).append(locus_tag)

        if not hits_by_query:
            log(f" ! No hits found for {group_name} (ID {group_id})", log_handle, also_print)

    return hits_by_query, best_hits_per_locus

######
# new section to run on folders of genomes
#####

# --- folder helper ---
def get_summary_info(gbk_file):
    record = next(SeqIO.parse(gbk_file, "genbank"))
    filename = Path(gbk_file).stem
    accession = record.id
    locus = record.name
    organism = record.annotations.get("organism", "-")
    definition = record.description
    return filename, accession, locus, organism, definition

# ---  Core Function ---
def run_wirehunter_on_file(gbk_file, query_file, args, out_dir) -> tuple:
    gbk_file = Path(gbk_file)
    prefix = gbk_file.stem
    temp_dir = tempfile.mkdtemp(dir=out_dir)

    output_gbk = out_dir / f"{prefix}_clusters.gbk"
    output_hits = out_dir / f"{prefix}_best_hits.tsv"
    log_path = out_dir / f"{prefix}_wirehunter.log"

    filename, accession, locus, organism, definition = get_summary_info(gbk_file)

    cluster_found = 0
    found_locus = "-"

    # Step 1: Parse query
    query_groups = read_query_tsv(query_file)

    # Step 2: Extract genome proteins
    proteins = extract_proteins_from_genbank(gbk_file)
    if not proteins:
        return accession, 0, locus, organism, definition

    prot_fasta = Path(temp_dir) / "genome_proteins.faa"
    write_fasta(proteins, prot_fasta)

    # Step 3: Build BLAST db
    subprocess.run([
        "makeblastdb", "-in", str(prot_fasta), "-dbtype", "prot",
        "-parse_seqids", "-out", str(Path(temp_dir) / "blastdb")
    ], check=True)

    # Map protein_id -> locus_tag
    locus_tag_map = {}
    for record in SeqIO.parse(gbk_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS" and "translation" in feature.qualifiers:
                locus_tag = feature.qualifiers.get("locus_tag", ["-"])[0]
                protein_id = feature.qualifiers.get("protein_id", [None])[0]
                if protein_id:
                    locus_tag_map[protein_id] = locus_tag

    # Step 4: Run BLAST for all query groups
    all_hits = {}
    all_best_hits_per_locus = {}
    with open(log_path, "w") as log_file:
        for group_id, entries in query_groups.items():
            group_fasta = Path(temp_dir) / f"query_{group_id}.faa"
            group_seqs = [
                SeqRecord(Seq(seq), id=f"{group_id}_{i}", description=name)
                for i, (name, seq, organism, protein_id) in enumerate(entries)
            ]
            write_fasta(group_seqs, group_fasta)

            result_xml = run_blastp(group_fasta, Path(temp_dir) / "blastdb", temp_dir)
            group_name = entries[0][0] if entries else f"group_{group_id}"

            hits, best_hits = parse_blast_results(
                result_xml, group_id, group_name, query_groups,
                args.min_coverage, log_file, True, args.min_identity, locus_tag_map
            )
            all_hits[group_id] = list(best_hits.keys())
            for locus, hit in best_hits.items():
                if locus not in all_best_hits_per_locus or hit["bit_score"] > all_best_hits_per_locus[locus]["bit_score"]:
                    all_best_hits_per_locus[locus] = hit

    # Step 5: Find cluster

    # Predefine default values in case no cluster is found
    total_genes = 0
    num_hits = 0
    hit_fraction = 0.0
    group_id_to_name = {
        group_id: entries[0][0] if entries else group_id
        for group_id, entries in query_groups.items()
    }
    missing_group_names = ", ".join(group_id_to_name[gid] for gid in sorted(query_groups.keys()))

    cluster_loci = find_flexible_cluster(
        hits_by_group=all_hits,
        min_hits=args.min_hits,
        gbk_path=gbk_file,
        max_distance=args.max_distance,
        window=args.window
    )

    if cluster_loci:
        cluster_found = 1
        extracted, region_locus_tags = extract_region_from_gbk(gbk_file, cluster_loci, args.window)
        found_locus = extracted[0].name if extracted else "-"

        with open(output_gbk, "w") as out:
            SeqIO.write(extracted, out, "genbank")

        total_genes = 0
        seen_query_names = set()  # track unique best_query_name values

        with open(output_hits, "w") as out:
            out.write("locus_tag\tbest_query_name\torganism\tprotein_id\tidentity\n")
            for locus in region_locus_tags:
                total_genes += 1
                best_hit = all_best_hits_per_locus.get(locus)
                if best_hit:
                    query_name = best_hit["query_name"]
                    seen_query_names.add(query_name)
                    out.write(f"{locus}\t{query_name}\t{best_hit['query_organism']}\t"
                              f"{best_hit['query_protein_id']}\t{best_hit['identity']:.1f}\n")
                else:
                    out.write(f"{locus}\tNo hit\t-\t-\t-\n")

        num_hits = len(seen_query_names)
        found_group_ids = set()
        for locus in region_locus_tags:
            best_hit = all_best_hits_per_locus.get(locus)
            if best_hit:
                group_id = best_hit["query_id"].split("_")[0]
                found_group_ids.add(group_id)

    hit_fraction = round(num_hits / total_genes, 2) if total_genes > 0 else 0.0

    all_group_ids = set(query_groups.keys())
    if cluster_found:
        missing_ids = all_group_ids - found_group_ids
        missing_names = [group_id_to_name.get(gid, gid) for gid in sorted(missing_ids)]
        missing_group_names = ", ".join(missing_names)
    else:
        missing_group_names = ""

    # Step 6: Find best protein ID where best_query_name is 'OmcE'
    
    best_interesting_protein_id = "-"
    if cluster_loci:
        for locus in cluster_loci:
            hit = all_best_hits_per_locus.get(locus)
            if hit and hit["query_name"] == "OmcE":
                best_interesting_protein_id = hit["query_protein_id"]
                break  # take first match
    
    # Delete temp folder unless user wants to keep it
    if not args.keep_temp and os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)

    return (filename, accession, cluster_found, found_locus, organism,
            definition, best_interesting_protein_id,
            total_genes, num_hits, hit_fraction, missing_group_names)

# --- New Main Entry Point ---
def main():
    args = parse_args()
    out_dir = Path(args.output)
    out_dir.mkdir(exist_ok=True)

    summary_path = out_dir / "cluster_summary.tsv"
    with open(summary_path, "w") as summary_out:
        summary_out.write("Filename\tAccession\tpresent\tLocus\tOrganism\tDefinition\tOmcE_protein_id\tnum_genes_in_region\tnum_hits\thit_fraction\tmissing_groups\n")
        
        if args.batch:
            gbk_files = list(Path(args.batch).rglob("*.gbk")) + list(Path(args.batch).rglob("*.gbff"))
            for gbk_file in gbk_files:
                print(f"--Processing {gbk_file.name}...")
                filename, accession, found, locus, org, desc, pore_id, total_genes, num_hits, hit_fraction, missing_groups = run_wirehunter_on_file(
                    gbk_file, args.query_file, args, out_dir
                )
                summary_out.write(f"{filename}\t{accession}\t{found}\t{locus}\t{org}\t{desc}\t{pore_id}\t"
                  f"{total_genes}\t{num_hits}\t{hit_fraction:.2f}\t{missing_groups}\n")

        elif args.genbank_file:
            filename, accession, found, locus, org, desc, pore_id, total_genes, num_hits, hit_fraction, missing_groups = run_wirehunter_on_file(
                args.genbank_file, args.query_file, args, out_dir
            )
            summary_out.write(f"{filename}\t{accession}\t{found}\t{locus}\t{org}\t{desc}\t{pore_id}\t"
                  f"{total_genes}\t{num_hits}\t{hit_fraction:.2f}\t{missing_groups}\n")

        else:
            print("❌ Error: You must provide either --genbank_file or --batch")
            sys.exit(1)

    print(f"✅ All done. Summary saved to {summary_path}\n")
    
if __name__ == "__main__":
    main()