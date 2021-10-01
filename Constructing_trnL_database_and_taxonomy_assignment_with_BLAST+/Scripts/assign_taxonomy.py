#!/usr/bin/env python3
"""Generates a CSV assigning taxonomy to sequences based on BLAST
   results."""

import sys

from Bio import SeqIO

IN_BLAST = sys.argv[1]
IN_HIT_FASTA = sys.argv[2]
IN_SEQ_FASTA = sys.argv[3]
IN_LEN_CUTOFF = float(sys.argv[4])
IN_CONS_CUTOFF = float(sys.argv[5])

def fasta_to_len_dict(in_fasta):
    """Returns a dictionary containing containing the length of each
       sequence as the value, with the sequence ID as the key."""
    len_dict = {}
    for rec in SeqIO.parse(in_fasta, 'fasta'):
        len_dict[rec.id] = len(rec.seq)
    return len_dict

def parse_blast_results(in_blast, in_len_cutoff, hit_len_dict, seq_len_dict):
    """Parses BLAST results for each sequence."""
    blast_dict = {}
    with open(in_blast, 'r') as in_handle:
        for line in in_handle:
            line = line.rstrip().split('\t')
            seq_id = line[0]
            seq_len = seq_len_dict[line[0]]
            if '{};'.format(line[1]) in hit_len_dict:
                hit_id = line[1].split('|')[0]
                hit_len = hit_len_dict['{};'.format(line[1])]
            else:
                hit_id = '{}.;'.format(line[1].split('|')[0])
                hit_len = hit_len_dict['{}.;'.format(line[1])]
            aln_idnt = float(line[2])
            aln_eval = -1*float(line[10])
            aln_min_frac = min((seq_len/hit_len), (hit_len/seq_len))
            if aln_min_frac >= in_len_cutoff:
                if seq_id not in blast_dict:
                    blast_dict[seq_id] = [set([hit_id]), aln_idnt, aln_eval,
                                          aln_min_frac]
                else:
                    if (
                        aln_idnt >= blast_dict[seq_id][1] and
                        aln_eval >= blast_dict[seq_id][2] and
                        aln_min_frac >= blast_dict[seq_id][3]
                    ):
                        if (
                            aln_idnt > blast_dict[seq_id][1] or
                            aln_eval > blast_dict[seq_id][2] or
                            aln_min_frac > blast_dict[seq_id][3]
                        ):
                            blast_dict[seq_id] = [set([hit_id]), aln_idnt,
                                                  aln_eval, aln_min_frac]
                        else:
                            blast_dict[seq_id][0].add(hit_id)
    return blast_dict

def find_cons(in_cons_cutoff, taxa):
    """Returns the best possible taxonomy assigment based on the cutoff
       value."""
    taxa = [t.split(';')[::-1] for t in taxa]
    tax_total = len(taxa)
    tax_cons = []
    min_cons = 0
    for tax_idx in range(len(taxa[0])):
        tax_count_dict = {}
        for tax in taxa:
            if tax[tax_idx] in tax_count_dict:
                tax_count_dict[tax[tax_idx]] += 1
            else:
                tax_count_dict[tax[tax_idx]] = 1
        max_tax = max(tax_count_dict, key=(lambda key: tax_count_dict[key]))
        max_count = tax_count_dict[max_tax]
        max_count_taxa = [t for t in tax_count_dict if tax_count_dict[t] ==
                          max_count]
        cons_val = max_count/tax_total
        if len(max_count_taxa) == 1 and cons_val >= in_cons_cutoff:
            tax_cons.append('"{}"'.format(max_tax))
            if not min_cons:
                min_cons = cons_val
        else:
            tax_cons.append('NA')
    return tax_cons[::-1], min_cons

def assign_tax(in_seq_fasta, in_cons_cutoff, blast_dict):
    """Assigns taxonomy based on the cutoff value and the BLAST
       results."""
    print(
        ','.join(
            [
                'Sequence', 'Kingdom', 'Phylum', 'Subphylum', 'Class', 'Order',
                'Family', 'Genus', 'Species', 'Consensus', 'Coverage',
                'Identity', 'E-value'
            ]
        )
    )
    for rec in SeqIO.parse(in_seq_fasta, 'fasta'):
        out_line = ['"{}"'.format(rec.seq)]
        if rec.id in blast_dict:
            taxa, aln_idnt, aln_eval, aln_min_frac = blast_dict[rec.id]
            tax, min_cons = find_cons(in_cons_cutoff, taxa)
            out_line.extend(tax)
            out_line.extend(str(v) for v in [min_cons, aln_min_frac, aln_idnt,
                                             -1*aln_eval])
        else:
            out_line.extend(['NA']*11)
        print(','.join(out_line))

HIT_LEN_DICT = fasta_to_len_dict(IN_HIT_FASTA)
SEQ_LEN_DICT = fasta_to_len_dict(IN_SEQ_FASTA)
BLAST_DICT = parse_blast_results(IN_BLAST, IN_LEN_CUTOFF, HIT_LEN_DICT,
                                 SEQ_LEN_DICT)
assign_tax(IN_SEQ_FASTA, IN_CONS_CUTOFF, BLAST_DICT)
