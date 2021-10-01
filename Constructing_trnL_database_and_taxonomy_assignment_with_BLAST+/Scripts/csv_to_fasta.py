#!/usr/bin/env python3
"""Converts an OTU table in the CSV format into a FASTA."""

import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

IN_CSV = sys.argv[1]

with open(IN_CSV, 'r') as in_handle:
    for line in in_handle:
        for seq_idx, seq in enumerate(line.rstrip().split(',')[1:]):
            rec = SeqRecord(Seq(seq), id=str(seq_idx+1), description='')
            SeqIO.write(rec, sys.stdout, 'fasta')
        break
