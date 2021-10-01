#!/usr/bin/env python3
"""Generates a FASTA file with taxonomy information from XML documents
   downloaded from NCBI."""

import sys
import xml.etree.ElementTree as ET

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

IN_SEQ_XML = sys.argv[1]
IN_REQ_RANK = sys.argv[2]
IN_EXC_NAME = sys.argv[3]
IN_TAX_XMLS = sys.argv[4:]

def seq_xml_to_dict(in_seq_xml):
    """Creates a sequence dictionary with the sequence as the value and
       the taxonomy ID as the key."""
    seq_dict = {}
    tree = ET.parse(in_seq_xml)
    root = tree.getroot()
    for tseq in root:
        taxid = tseq.find('TSeq_taxid').text
        accver = tseq.find('TSeq_accver').text
        seq = tseq.find('TSeq_sequence').text
        if taxid not in seq_dict:
            seq_dict[taxid] = [(accver, seq)]
        else:
            seq_dict[taxid].append((accver, seq))
    return seq_dict

def format_name(sciname):
    """Reformats a scientific name to remove disallowed characters."""
    return sciname.replace(' ', '_').replace(';', ',').replace('|', '_')

def tax_xmls_to_dict(in_tax_xmls, in_req_rank, in_exc_name):
    """Loads taxonomy information into a dictionary with the taxonomy ID
       as the key and the taxonomy as the value."""
    acc_ranks = ['kingdom', 'phylum', 'subphylum', 'class', 'order', 'family',
                 'genus', 'species']
    tax_dict = {}
    for in_tax_xml in in_tax_xmls:
        tree = ET.parse(in_tax_xml)
        root = tree.getroot()
        for taxon in root:
            taxon_dict = {}
            taxid = taxon.find('TaxId').text
            rank = taxon.find('Rank').text
            sciname = format_name(taxon.find('ScientificName').text)
            if rank in acc_ranks:
                taxon_dict[rank] = sciname
            for taxon_ex in taxon.find('LineageEx'):
                rank_ex = taxon_ex.find('Rank').text
                sciname_ex = format_name(taxon_ex.find('ScientificName').text)
                if rank_ex in acc_ranks:
                    taxon_dict[rank_ex] = sciname_ex
            if (in_req_rank in taxon_dict and
                not any(taxon_dict[r] == in_exc_name for r in taxon_dict)):
                ranks = []
                max_idx = max(acc_ranks.index(a) for a in acc_ranks)
                for acc_idx, acc_rank in enumerate(acc_ranks):
                    if acc_rank in taxon_dict:
                        ranks.append(taxon_dict[acc_rank])
                    elif acc_idx < max_idx:
                        ranks.append('Unclassified')
                tax_dict[taxid] = ranks
                aka_taxids = taxon.find('AkaTaxIds')
                if aka_taxids:
                    for aka_taxid in aka_taxids:
                        tax_dict[aka_taxid.text] = ranks
    return tax_dict

def dicts_to_fasta(seq_dict, tax_dict):
    """Generates a FASTA from information contained in the sequence and
       taxonomy dictionaries."""
    for taxid in sorted(tax_dict, key=lambda k: tax_dict[k]):
        if taxid in seq_dict:
            for accver, seq in seq_dict[taxid]:
                rec = SeqRecord(Seq(seq), id=';'.join(tax_dict[taxid]),
                                description=accver)
                SeqIO.write(rec, sys.stdout, 'fasta')
                rec.seq = rec.seq.reverse_complement()
                SeqIO.write(rec, sys.stdout, 'fasta')

SEQ_DICT = seq_xml_to_dict(IN_SEQ_XML)
TAX_DICT = tax_xmls_to_dict(IN_TAX_XMLS, IN_REQ_RANK, IN_EXC_NAME)
dicts_to_fasta(SEQ_DICT, TAX_DICT)
