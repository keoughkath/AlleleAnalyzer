#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pam_pos_genome.py finds all PAM positions for given PAMs in a given genome.
Written in Python v 3.6.1.
Kathleen Keough et al 2018.

Usage:
    pam_pos_genome.py <chrom> <fasta> <cas_list> <out>

Arguments:
    chrom             Chromosome being analyzed.
    fasta             Fasta file for genome being analyzed.
    cas_list          Comma separated (no spaces!) list of Cas varieties to evaluate, options below.
    out               Out prefix for returned files.

Available Cas types = cpf1,SpCas9,SpCas9_VRER,SpCas9_EQR,SpCas9_VQR_1,SpCas9_VQR_2,StCas9,StCas9_2,SaCas9,SaCas9_KKH,nmCas9,cjCas9
"""

import numpy as np
import sys, os
from pyfaidx import Fasta
import pandas as pd
import regex
import re
from Bio import SeqIO
from docopt import docopt

__version__='1.0.0'

# 3 and 5 prime cas lists

# Get absolute path for gen_targ_dfs.py, and edit it for cas_object.py
ef_path = os.path.dirname(os.path.realpath(__file__))
cas_obj_path = ef_path.replace('preprocessing/find_pams_in_reference','scripts/')
sys.path.append(cas_obj_path)

metadata_path = ef_path.replace('find_pams_in_reference','')
sys.path.append(metadata_path)
# Import cas_object
import cas_object as cas_obj
from get_metadata import add_metadata

# get rid of annoying false positive Pandas error

pd.options.mode.chained_assignment = None

# three-prime PAMs, e.g. Cas9, PAM is 3' of the sgRNA sequence
tpp_for = {}

tpp_for['SpCas9'] = r'[atcg]gg' # SpCas9, SpCas9-HF1, eSpCas1.1
tpp_for['SpCas9_VRER'] = r'[atcg]gcg' # SpCas9 VRER variant
tpp_for['SpCas9_EQR'] = r'[actg]gag' # SpCas9 EQR variant
tpp_for['SpCas9_VQR_1'] = r'[atcg]ga' # SpCas9 VQR variant 1
tpp_for['SpCas9_VQR_2'] = r'[atcg]g[atcg]g' # SpCas9 VQR variant 2
tpp_for['StCas9'] = r'[actg]{2}agaa' # S. thermophilus Cas9
tpp_for['StCas9_2'] = r'[actg]gg[actg]g' # S. thermophilus Cas9 2
tpp_for['SaCas9'] = r'[atcg]{2}g[ag]{2}t' # SaCas9
tpp_for['SaCas9_KKH'] = r'[atcg]{3}[ag]{2}t' # SaCas9 KKH variant
tpp_for['nmCas9'] = r'[atcg]{4}g[ac]tt' # nmCas9
tpp_for['cjCas9'] = r'[actg]{4}aca' # campylobacter jejuni Cas9

# find 3' PAMs on antisense strand (reverse complement)
tpp_rev = {}

tpp_rev['SpCas9_rev'] = r'cc[atcg]' # SpCas9 reverse complement 
tpp_rev['SpCas9_VRER_rev'] = r'cgc[atcg]' # SpCas9 VRER variant reverse complement 
tpp_rev['SpCas9_EQR_rev'] = r'ctc[actg]' # SpCas9 EQR variant reverse complement 
tpp_rev['SpCas9_VQR_1_rev'] = r'[atcg]tc[atcg]' # SpCas9 VQR variant 1 reverse complement 
tpp_rev['SpCas9_VQR_2_rev'] = r'c[atcg]c[atcg]' # SpCas9 VQR variant 2 reverse complement 
tpp_rev['StCas9_rev'] = r'ttct[actg]{2}' # S. thermophilus Cas9 reverse complement
tpp_rev['StCas9_2_rev'] = r'g[atcg]gg[atcg]' # S. thermophilus Cas9 2 reverse complement 
tpp_rev['SaCas9_rev'] = r'a[tc]{2}c[atcg]{2}' # SaCas9 reverse complement 
tpp_rev['SaCas9_KKH_rev'] = r'a[tc]{2}[atcg]{3}' # SaCas9 KKH variant reverse complement 
tpp_rev['nmCas9_rev'] = r'aa[tg]c[atcg]{4}' # NmCas9 reverse complement 
tpp_rev['cjCas9_rev'] = r'tgt[actg]{4}' # campylobacter jejuni Cas9

# five-prime PAMs, e.g. cpf1, PAM is 5' of the sgRNA sequence
fpp_for = {}

fpp_for['cpf1'] = r'ttt[atcg]' # Cpf1, PAM 5' of guide

# find 5' PAMs on antisense strand (reverse complement)
fpp_rev = {}

fpp_rev['cpf1_rev'] = r'[atcg]aaa' # Cpf1, PAM 5' of guide

def find_the_pams(seqio_object):

    # this will be the ultimate output, a dictionary linking PAM names to their positions

    pam_positions = {}

    # get sequence of inputted seqIO object

    sequence = str(seqio_object.seq)
    seqid = seqio_object.id

    # set inputs to variables, string for sequence of region of interest
    # ends up being region_seq. 

    chrom = seqid.split(':')[0]
    low = int(seqid.split(':')[1].split('-')[0]) # low end of region of interest
    high = int(seqid.split(':')[1].split('-')[1]) # high end of region of interest

    # interior function to get PAM sites for each PAM motif
    
    def get_pam_starts(pam_regex,sequence):
        starts = set()
        for pam in regex.finditer(pam_regex, sequence,regex.IGNORECASE,overlapped=True):
            starts.add(pam.end()+low+1) 
        return(set(starts))

    for key,value in pam_dict.items():
        pam_positions[key] = get_pam_starts(value,sequence)

    # return dictionary of PAMs and their positions in the sequence

    return(pam_positions)

def find_spec_pams(cas_obj,python_string,orient="3'"):
    # orient specifies whether this is a 3prime PAM (e.g. Cas9, PAM seq 3' of sgRNA)
    # or a 5prime PAM (e.g. cpf1, PAM 5' of sgRNA)

    # get sequence 

    sequence = python_string

    # get PAM sites (the five prime three prime thing will need to be reversed for cpf1)

    def get_pam_fiveprime(pam_regex,sequence):
        starts = set()
        for pam in regex.finditer(pam_regex, sequence,regex.IGNORECASE,overlapped=True):
            starts.add(pam.start()+1) 
        return(set(starts))

    def get_pam_threeprime(pam_regex,sequence):
        starts = set()
        for pam in regex.finditer(pam_regex, sequence,regex.IGNORECASE,overlapped=True):
            starts.add(pam.end()) 
        return(set(starts))

    if orient == "3'":
        for_starts = get_pam_fiveprime(cas_obj.forwardPam_regex() ,sequence)
        rev_starts = get_pam_threeprime(cas_obj.reversePam_regex(),sequence)
    elif orient == "5'":
        for_starts = get_pam_threeprime(cas_obj.forwardPam_regex(),sequence)
        rev_starts = get_pam_fiveprime(cas_obj.reversePam_regex(),sequence)

    return(for_starts,rev_starts)

def main(args):

    # keep track of chromosome since this will be run with bash script

    chrom = args['<chrom>']

    # define output prefix

    outprefix = args['<out>']

    # make Fasta object for genome of choice, e.g. hg19

    genome = Fasta(args['<fasta>'],as_raw=True)

    cas_list = args['<cas_list>'].split(',')

    # get set of positions for each type of cas

    for cas in cas_list:
        current_cas = cas_obj.get_cas_enzyme(cas, os.path.join(cas_obj_path,'CAS_LIST.txt'))
        for_starts, rev_starts = find_spec_pams(current_cas,str(genome[str(chrom)]), orient=current_cas.primeness)
        savestr_for = f'{outprefix}'+str(chrom)+'_'+str(cas) + '_pam_sites_for'
        savestr_rev = f'{outprefix}'+str(chrom)+'_'+str(cas) + '_pam_sites_rev'
        np.save(savestr_for + '.npy',list(for_starts))
        np.save(savestr_rev + '.npy',list(rev_starts))
        np.savetxt(savestr_for + '.txt',list(for_starts))
        np.savetxt(savestr_rev + '.txt',list(rev_starts))


if __name__ == '__main__':
    arguments = docopt(__doc__, version=__version__)
    main(arguments)
