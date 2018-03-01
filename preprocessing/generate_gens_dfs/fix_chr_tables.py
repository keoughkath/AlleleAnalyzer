#!/usr/bin/env python
#-*- coding: utf-8 -*-

"""
bcftools doesn't split every column, so this finishes the job!
Intended to be used with get_chr_tables.sh.
Written in Python version 3.6.1.
Kathleen Keough 2017.
"""

import pandas as pd
import sys
import os
import regex as re

__version__ = '0.0.0'


def fix_multiallelics(cell):
    """
    bcftools doesn't complete splitting of multiallelic variant sites.
    :param cell: genotype, str.
    :return: genotype as is if not multiallelic otherwise split multiallelic genotype, str.
    """
    splitters = [',', ';']
    if any(splitter in str(cell) for splitter in splitters):
        cell = re.split(';|,', cell)[0]
    return cell


def main():
    vars = pd.read_csv(sys.argv[1], sep='\t', header=None, names=['chrom', 'pos', 'ref', 'alt', 'genotype'],
                       usecols=['chrom', 'pos', 'ref', 'alt', 'genotype'])
    chrom = sys.argv[2]
    if 'chr' in str(vars.chrom.iloc[0]):
        vars['chrom'] = vars['chrom'].str[3:]
    vars = vars.query('chrom == @chrom')
    vars_fixed = vars.applymap(fix_multiallelics)
    
    vars_fixed.to_hdf(os.path.join(sys.argv[3], f'chr{chrom}_gens.hdf5'), 'all', complib='blosc')


if __name__ == '__main__':
    main()
