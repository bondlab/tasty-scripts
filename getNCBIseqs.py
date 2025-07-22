#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 28 11:51:05 2025

@author: daniel
"""

from Bio import Entrez, SeqIO

# Use your email here as required by NCBI
Entrez.email = "dbond@umn.edu"

ids = [
    "WP_161598787.1", "WP_004512831.1", "WP_214186954.1", "WP_039743990.1",
    "WP_241426340.1", "WP_214175778.1", "WP_215733991.1", "WP_237560669.1",
    "WP_011940632.1", "WP_012646522.1", "WP_214299667.1", "WP_243689285.1",
    "WP_243689284.1", "WP_214170156.1", "WP_145017695.1"
]

with Entrez.efetch(db="protein", id=",".join(ids), rettype="fasta", retmode="text") as handle:
    records = list(SeqIO.parse(handle, "fasta"))

with open("OMcEBB.fasta", "w") as out_f:
    SeqIO.write(records, out_f, "fasta")