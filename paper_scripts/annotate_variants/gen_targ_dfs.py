#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
gen_targ_dfs.py generates necessary variant targetability data for use with ExcisionFinder. Written in Python v 3.6.1.
Kathleen Keough 2017.

Usage:
    gen_targ_dfs.py <gens_file> <chrom> <pams_dir> <hg19_fasta> <out_dir>

Arguments:
    gens_file           Explicit genotypes file generated by get_chr_tables.sh
    chrom               Chromosome in the file.
    pams_dir            Directory where pam locations in hg19 are located. These are pregenerated and on the EF Github.
    hg19_fasta          Fasta file for hg19 (GrCh37.fa provided on EF github)
    out_dir             Directory in which to save the output files.
"""

import pandas as pd
import numpy as np
from docopt import docopt
import os
import crisprtools
from pyfaidx import Fasta

__version__ = '0.0.0'

# 3 and 5 prime cas lists

TP_CAS_LIST = ['SpCas9', 'SpCas9_VRER', 'SpCas9_EQR', 'SpCas9_VQR_1',
               'SpCas9_VQR_2', 'StCas9', 'StCas9_2', 'SaCas9', 'SaCas9_KKH', 'nmCas9', 'cjCas9']

FP_CAS_LIST = ['cpf1']

CAS_LIST = ['cpf1', 'SpCas9', 'SpCas9_VRER', 'SpCas9_EQR', 'SpCas9_VQR_1',
            'SpCas9_VQR_2', 'StCas9', 'StCas9_2', 'SaCas9', 'SaCas9_KKH', 'nmCas9', 'cjCas9']


def get_range_upstream(pam_pos):
    """
    Get positions 20 bp upstream, i.e. for forward 3' PAMs or reverse 5' PAMs
    :param pam_pos: position of PAM, int.
    :return: sgRNA seed region positions, set of ints.
    """
    sgrna = set(range(pam_pos - 21, pam_pos))
    return sgrna


def get_range_downstream(pam_pos):
    """
    Get positions 20 bp upstream, i.e. for forward 3' PAMs or reverse 5' PAMs
    :param pam_pos: position of PAM, int.
    :return: sgRNA seed region positions, set of ints.
    """
    sgrna = set(range(pam_pos + 1, pam_pos + 21))
    return sgrna


def makes_breaks_pam(chrom, pos, ref, alt, hg19):
    """
    Determine if cas in question makes or breaks PAM sites.
    :param chrom: chromosome, int.
    :param pos: position, int.
    :param ref: ref genotype, str.
    :param alt: alt genotype, str.
    :param hg19: hg19 fasta file, fasta.
    :return:
    """
    makes_pam = [0] * len(CAS_LIST)
    breaks_pam = [0] * len(CAS_LIST)

    if '<' in alt:
        return makes_pam + breaks_pam

    # if alt is not a special case (CNV or SV), continue checking the new sequence

    ref_seq = hg19[str(chrom)][pos - 11:pos + 10]

    if len(ref) > len(alt):  # handles deletions
        alt_seq = hg19[str(chrom)][pos - 11:pos - 1] + alt + hg19[str(chrom)][
                                                             pos + len(ref) + len(alt) - 3:pos + len(ref) + len(
                                                                 alt) - 3 + 10]
    else:
        alt_seq = hg19[str(chrom)][pos - 11:pos - 1] + alt + hg19[str(chrom)][
                                                             pos + len(alt) - 2:pos + len(alt) - 2 + 10]

    counter = -1
    for cas in CAS_LIST:
        # print(cas,ref,alt,pos)
        counter += 1
        if cas == 'cpf1':
            ref_pams_for, ref_pams_rev = crisprtools.find_spec_pams(cas, ref_seq, orient='5prime')
            alt_pams_for, alt_pams_rev = crisprtools.find_spec_pams(cas, alt_seq, orient='5prime')
        else:
            ref_pams_for, ref_pams_rev = crisprtools.find_spec_pams(cas, ref_seq)
            alt_pams_for, alt_pams_rev = crisprtools.find_spec_pams(cas, alt_seq)

        if len(alt_pams_for) - len(ref_pams_for) > 0 or len(alt_pams_rev) - len(ref_pams_rev) > 0:
            makes_pam[counter] = 1
        elif len(ref_pams_for) - len(alt_pams_for) > 0 or len(ref_pams_rev) - len(alt_pams_rev) > 0:
            breaks_pam[counter] = 1
    return makes_pam + breaks_pam


def get_made_broke_pams(df, chrom, hg19):
    """
    Apply makes_breaks_pams to a df.
    :param df: gens df generated by get_chr_tables.sh, available on EF github.
    :param chrom: chromosome currently being analyzed.
    :param hg19: hg19 fasta, pyfaidx format.
    :return: dataframe with indicators for whether each variant makes/breaks PAMs, pd df.
    """
    df['makes_cpf1'], df['makes_SpCas9'], df['makes_SpCas9_VRER'], df['makes_SpCas9_EQR'], \
        df['makes_SpCas9_VQR_1'], df['makes_SpCas9_VQR_2'], df['makes_StCas9'], df['makes_StCas9_2'], \
        df['makes_SaCas9'], df['makes_SaCas9_KKH'], df['makes_nmCas9'], df['makes_cjCas9'], \
        df['breaks_cpf1'], df['breaks_SpCas9'], df['breaks_SpCas9_VRER'], df['breaks_SpCas9_EQR'], \
        df['breaks_SpCas9_VQR_1'], df['breaks_SpCas9_VQR_2'], df['breaks_StCas9'], df['breaks_StCas9_2'], \
        df['breaks_SaCas9'], df['breaks_SaCas9_KKH'], df['breaks_nmCas9'], df['breaks_cjCas9'] = \
        zip(*df.apply(lambda row: makes_breaks_pam(chrom, row['pos'], row['ref'], row['alt'], hg19), axis=1))
    return df


def main(args):
    out_dir = args['<out_dir>']
    pams_dir = args['<pams_dir>']
    gens = args['<gens_file>']
    chrom = args['<chrom>']
    hg19 = Fasta(args['<hg19_fasta>'], as_raw=True)

    os.makedirs(out_dir, exist_ok=True)

    # load chromosome variants

    gens = pd.read_hdf(gens, 'all')
    chr_variants = set(gens['pos'].tolist())

    # save locations of PAM proximal variants to dictionary

    pam_prox_vars = {}

    # get variants within sgRNA region for 3 prime PAMs (20 bp upstream of for pos and vice versa)

    for cas in TP_CAS_LIST:
        print(cas)
        cas_prox_vars = []
        pam_dict = {}
        pam_for_pos = np.load(os.path.join(pams_dir, f'chr{chrom}_{cas}_pam_sites_for.npy')).tolist()
        pam_rev_pos = np.load(os.path.join(pams_dir, f'chr{chrom}_{cas}_pam_sites_rev.npy')).tolist()
        for pos in pam_for_pos:
            prox_vars = set(get_range_upstream(pos)) & chr_variants
            cas_prox_vars.extend(prox_vars)
            pam_dict[pos] = prox_vars
        for pos in pam_rev_pos:
            prox_vars = set(get_range_downstream(pos)) & chr_variants
            cas_prox_vars.extend(prox_vars)
            pam_dict[pos] = prox_vars
        pam_prox_vars[cas] = cas_prox_vars

    # same for five prime pams

    for cas in FP_CAS_LIST:
        print(cas)
        cas_prox_vars = []
        pam_dict = {}
        pam_for_pos = np.load(os.path.join(pams_dir, f'chr{chrom}_{cas}_pam_sites_for.npy')).tolist()
        pam_rev_pos = np.load(os.path.join(pams_dir, f'chr{chrom}_{cas}_pam_sites_rev.npy')).tolist()
        for pos in pam_for_pos:
            prox_vars = set(get_range_downstream(pos)) & chr_variants
            cas_prox_vars.extend(prox_vars)
            pam_dict[pos] = prox_vars
        for pos in pam_rev_pos:
            prox_vars = set(get_range_upstream(pos)) & chr_variants
            cas_prox_vars.extend(prox_vars)
            pam_dict[pos] = prox_vars
        pam_prox_vars[cas] = cas_prox_vars

    chrdf = get_made_broke_pams(gens, chrom, hg19)
    # make_break_df.to_hdf(os.path.join(out_dir, f'chr{chrom}_make_break.hdf5'), 'all', complib='blosc')

    for cas in CAS_LIST:
        print(cas)
        spec_pam_prox_vars = pam_prox_vars[cas]
        chrdf[f'var_near_{cas}'] = chrdf['pos'].isin(spec_pam_prox_vars).astype(int)

    chrdf.to_hdf(os.path.join(out_dir, f'chr{chrom}_targ.hdf5'), 'all', mode='w', format='table', data_columns=True,
                 chunksize=200000, complib='blosc')


if __name__ == '__main__':
    arguments = docopt(__doc__, version=__version__)
    main(arguments)
