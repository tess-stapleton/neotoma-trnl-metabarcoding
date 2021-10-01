#!/usr/bin/env python3
"""Generates a database from a Cutadapt-trimmed FASTA."""

import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

IN_FASTA = sys.argv[1]
IN_PRIMER_1 = sys.argv[2]
IN_PRIMER_2 = sys.argv[3]
MAX_LEN = int(sys.argv[4])
MIN_LEN = int(sys.argv[5])

def trim_seq(seq, in_primer_1, in_primer_2, max_len, min_len):
    """Returns the trimmed portion of a sequence provided it contains no
       ambiguous nucleotides, adapters appear to be complete, and the
       sequence of interest is in the correct size range."""
    acc_nucs = {'A', 'C', 'G', 'T'}
    primer_1_len = len(in_primer_1)
    primer_2_len = len(in_primer_2)
    trimmed_seq = ''
    ranges = [0, 0, 0]
    for nuc in seq:
        if nuc.islower(): 
            if ranges[1] == 0:
                ranges[0] += 1
            else:
                ranges[2] += 1
        else:
            ranges[1] += 1
    if (ranges[0] >= primer_1_len and ranges[1] <= max_len and
        ranges[1] >= min_len and ranges[2] >= primer_2_len):
        seq_start = ranges[0]-primer_1_len
        seq_end = ranges[0]+ranges[1]+primer_2_len
        seq_primers = seq[seq_start:seq_end]
        if set(seq_primers.upper()).issubset(acc_nucs):
            trimmed_seq = seq_primers[primer_1_len:-primer_2_len]
    return trimmed_seq

def cutadapt_fasta_to_dict(in_fasta, in_primer_1, in_primer_2, max_len,
                           min_len):
    """Builds a dictionary of trimmed sequences from the FASTA file with
       the taxonomy set as the key and sequences stored as values in
       sets."""
    fasta_db = {}
    for rec in SeqIO.parse(in_fasta, 'fasta'):
        trimmed_seq = trim_seq(str(rec.seq), in_primer_1, in_primer_2, max_len,
                               min_len)
        if trimmed_seq:
            tax = rec.id.split(' ')[0]
            if tax not in fasta_db:
                fasta_db[tax] = set()
            fasta_db[tax].add(trimmed_seq)
    return fasta_db

def dict_to_db_fasta(fasta_db):
    """Builds a database FASTA from the trimmed sequence dictionary."""
    for sorted_tax in sorted(t.split(';') for t in fasta_db):
        tax = ';'.join(sorted_tax)
        suffix = True if len(fasta_db[tax]) > 1 else False
        for seq_idx, seq in enumerate(sorted(fasta_db[tax])):
            if suffix:
                rec_id = '{}|{};'.format(tax, seq_idx+1)
            else:
                rec_id = '{};'.format(tax)
            rec = SeqRecord(Seq(seq), id=rec_id, description = '')
            SeqIO.write(rec, sys.stdout, 'fasta')

FASTA_DB = cutadapt_fasta_to_dict(IN_FASTA, IN_PRIMER_1, IN_PRIMER_2, MAX_LEN,
                                  MIN_LEN)
dict_to_db_fasta(FASTA_DB)
