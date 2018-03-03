# import modules

import pandas as pd
import regex
import re
from Bio import SeqIO

# get rid of annoying false positive Pandas error

pd.options.mode.chained_assignment = None

# "three-prime PAMs, e.g. Cas9, PAM is 3' of the sgRNA sequence"

# tuple which includes length intended to make it easier to generate 
# IGV files later than also show PAM site, not to be confused with 
# the manuscript analysis that uses only the non-N bases of PAMs 
# and analyses how that influences occurrence in the genome

tpp_for = {}

tpp_for['SpCas9'] = (r'[atcg]gg', 3) # SpCas9, SpCas9-HF1, eSpCas1.1
tpp_for['SpCas9_VRER'] = (r'[atcg]gcg', 4) # SpCas9 VRER variant
tpp_for['SpCas9_EQR'] = (r'[actg]gag', 4) # SpCas9 EQR variant
tpp_for['SpCas9_VQR_1'] = (r'[atcg]ga[atcg]', 4) # SpCas9 VQR variant 1
tpp_for['SpCas9_VQR_2'] = (r'[atcg]g[atcg]g', 4) # SpCas9 VQR variant 2
tpp_for['StCas9'] = (r'[actg]{2}agaa', 6) # S. thermophilus Cas9
tpp_for['StCas9_2'] = (r'[actg]gg[actg]g', 5) # S. thermophilus Cas9 2
tpp_for['SaCas9'] = (r'[atcg]{2}g[ag]{2}t', 6) # SaCas9
tpp_for['SaCas9_KKH'] = (r'[atcg]{3}[ag]{2}t', 6) # SaCas9 KKH variant
tpp_for['nmCas9'] = (r'[atcg]{4}g[ac]tt', 8) # nmCas9
tpp_for['cjCas9'] = (r'[actg]{4}aca', 7) # campylobacter jejuni Cas9

# find 3' PAMs on antisense strand (reverse complement)
tpp_rev = {}

tpp_rev['SpCas9_rev'] = (r'cc[atcg]', 3) # SpCas9 reverse complement 
tpp_rev['SpCas9_VRER_rev'] = (r'cgc[atcg]', 4) # SpCas9 VRER variant reverse complement 
tpp_rev['SpCas9_EQR_rev'] = (r'ctc[actg]', 4) # SpCas9 EQR variant reverse complement 
tpp_rev['SpCas9_VQR_1_rev'] = (r'[atcg]tc[atcg]', 4) # SpCas9 VQR variant 1 reverse complement <--
tpp_rev['SpCas9_VQR_2_rev'] = (r'c[atcg]c[atcg]', 4) # SpCas9 VQR variant 2 reverse complement 
tpp_rev['StCas9_rev'] = (r'ttct[actg]{2}', 6) # S. thermophilus Cas9 reverse complement
tpp_rev['StCas9_2_rev'] = (r'g[atcg]gg[atcg]', 5) # S. thermophilus Cas9 2 reverse complement 
tpp_rev['SaCas9_rev'] = (r't[tc]{2}c[atcg]{2}', 6) # SaCas9 reverse complement 
tpp_rev['SaCas9_KKH_rev'] = (r'a[tc]{2}[atcg]{3}', 6) # SaCas9 KKH variant reverse complement 
tpp_rev['nmCas9_rev'] = (r'aa[tg]c[atcg]{4}', 8) # NmCas9 reverse complement 
tpp_rev['cjCas9_rev'] = (r'tgt[actg]{4}', 7) # campylobacter jejuni Cas9
# "five-prime PAMs, e.g. cpf1, PAM is 5' of the sgRNA sequence"
fpp_for = {}

fpp_for['cpf1'] = (r'ttt[atcg]', 4) # Cpf1, PAM 5' of guide

# find 5' PAMs on antisense strand (reverse complement)
fpp_rev = {}

fpp_rev['cpf1_rev'] = (r'[atcg]aaa', 4) # Cpf1, PAM 5' of guide

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

def find_spec_pams(cas,python_string,orient='3prime'):
	# orient specifies whether this is a 3prime PAM (e.g. Cas9, PAM seq 3' of sgRNA)
	# or a 5prime PAM (e.g. cpf1, PAM 5' of sgRNA)

    # get sequence 

    sequence = python_string

    # get PAM sites (the five prime three prime thing will need to be reversed for cpf1)

    def get_pam_fiveprime(pam_regex,sequence):
        starts = []
        for pam in regex.finditer(pam_regex, sequence,regex.IGNORECASE,overlapped=True):
            starts.append(pam.start()) 
        return(starts)

    def get_pam_threeprime(pam_regex,sequence):
        starts = []
        for pam in regex.finditer(pam_regex, sequence,regex.IGNORECASE,overlapped=True):
            starts.append(pam.end()) 
        return(starts)

    if orient == '3prime':
        for_starts = get_pam_fiveprime(tpp_for[cas][0],sequence)
        rev_starts = get_pam_threeprime(tpp_rev[cas+'_rev'][0],sequence)
    elif orient == '5prime':
        for_starts = get_pam_threeprime(fpp_for[cas][0],sequence)
        rev_starts = get_pam_fiveprime(fpp_rev[cas+'_rev'][0],sequence)

    return(for_starts,rev_starts)