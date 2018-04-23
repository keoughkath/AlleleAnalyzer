#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script generates an arcplot demonstrating shared allele-specific targetable sites in a particular gene.
Kathleen Keough & Svetlana Lyalina 2017-2018.

Usage:
        gen_arcplot_df.py --type=<type> <gene> <filt> <annots> <cas> <hap_targs> <out_dir>

Arguments:
    type            type of df outputted, a for arcplot and s for set cover
    gene            gene
    filt            minimum % (integer) of included individuals that must share the outputted site pairs (0 for set cover)
    annots          gene annotations file
    cas             cas variety
    hap_targs       outputted by ExcisionFinder
    out_dir         prefix/directory to save the file
"""
import pandas as pd
import itertools
import ExcisionFinder as ef
from docopt import docopt
import numpy as np


def ppl_w_pair(df, pos1, pos2, cas):
    query1 = f'pos == @pos1 and hap1_{cas}'
    ppl_targ1_hap1 = df.query(query1)['sample'].tolist()
    query2 = f'pos == @pos2 and hap1_{cas}'
    ppl_targ2_hap1 = df.query(query2)['sample'].tolist()
    query3 = f'pos == @pos1 and hap2_{cas}'
    ppl_targ1_hap2 = df.query(query3)['sample'].tolist()
    query4 = f'pos == @pos2 and hap2_{cas}'
    ppl_targ2_hap2 = df.query(query4)['sample'].tolist()
    ppl_targ_hap1 = set(ppl_targ1_hap1).intersection(set(ppl_targ2_hap1))
    ppl_targ_hap2 = set(ppl_targ1_hap2).intersection(set(ppl_targ2_hap2))
    ppl_targ = ppl_targ_hap1.union(ppl_targ_hap2)
    return len(ppl_targ)


def ppl_ids_w_pair(df, pos1, pos2, cas):
    query1 = f'pos == @pos1 and hap1_{cas}'
    ppl_targ1_hap1 = df.query(query1)['sample'].tolist()
    query2 = f'pos == @pos2 and hap1_{cas}'
    ppl_targ2_hap1 = df.query(query2)['sample'].tolist()
    query3 = f'pos == @pos1 and hap2_{cas}'
    ppl_targ1_hap2 = df.query(query3)['sample'].tolist()
    query4 = f'pos == @pos2 and hap2_{cas}'
    ppl_targ2_hap2 = df.query(query4)['sample'].tolist()
    ppl_targ_hap1 = set(ppl_targ1_hap1).intersection(set(ppl_targ2_hap1))
    ppl_targ_hap2 = set(ppl_targ1_hap2).intersection(set(ppl_targ2_hap2))
    ppl_targ = ppl_targ_hap1.union(ppl_targ_hap2)
    return ppl_targ


def get_arcplot_df(indf, geneid, gene, cas, annots, out_dir):
    exons, n_coding, genestart, genestop, coding_positions, coding_exon_starts \
        = ef.get_coding_exons(str(geneid), annots, gene, out_dir)
    if cas == 'all':
        hap1_cols = list(filter(lambda x: x.startswith('hap1'), indf.columns))
        hap2_cols = list(filter(lambda x: x.startswith('hap2'), indf.columns))
        indf[f'hap1_all'] = indf[hap1_cols].any(axis=1)
        indf[f'hap2_all'] = indf[hap2_cols].any(axis=1)
    query = f'hap1_{cas} or hap2_{cas}'
    n_inds = len(indf['sample'].drop_duplicates().tolist())
    poisdf = indf.query(query)[['pos', 'sample', f'hap1_{cas}', f'hap2_{cas}']]  # positions of interest
    variants = sorted(poisdf.pos.drop_duplicates().tolist())
    variant1 = []
    variant2 = []
    for var1, var2 in itertools.combinations(variants, r=2):
        if (var1 != var2) and (max([var1, var2]) <= min([var1, var2]) + 10000) and (ef.targ_pair(var1, var2,
                                                                                                 coding_positions,
                                                                                                 coding_exon_starts)):
            variant1.append(var1)
            variant2.append(var2)
        else:
            continue

    targ_pairs_df = pd.DataFrame({'var1': variant1, 'var2': variant2})
    targ_pairs_df['n_inds'] = targ_pairs_df.apply(lambda row: ppl_w_pair(indf, row['var1'], row['var2'], cas), axis=1)
    targ_pairs_df.query('n_inds > 0', inplace=True)
    targ_pairs_df['percent_pop_covered'] = (targ_pairs_df['n_inds'] / n_inds) * 100.0

    return targ_pairs_df


def get_set_cover_df(indf, geneid, gene, cas, annots, out_dir):
    exons, n_coding, genestart, genestop, coding_positions, coding_exon_starts \
        = ef.get_coding_exons(str(geneid), annots, gene, out_dir)
    if cas == 'all':
        hap1_cols = list(filter(lambda x: x.startswith('hap1'), indf.columns))
        hap2_cols = list(filter(lambda x: x.startswith('hap2'), indf.columns))
        indf[f'hap1_all'] = indf[hap1_cols].any(axis=1)
        indf[f'hap2_all'] = indf[hap2_cols].any(axis=1)
    query = f'hap1_{cas} or hap2_{cas}'
    poisdf = indf.query(query)[['pos', 'sample', f'hap1_{cas}', f'hap2_{cas}']]  # positions of interest
    variants = sorted(poisdf.pos.drop_duplicates().tolist())
    variant1 = []
    variant2 = []
    inds = indf['sample'].drop_duplicates()
    # inds = pd.Series(np.load('all_samples.npy').tolist())
    ppl_included_dict = {}
    counter = -1
    for var1, var2 in itertools.combinations(variants, r=2):
        if (var1 != var2) and (max([var1, var2]) <= min([var1, var2]) + 10000) and (ef.targ_pair(var1, var2,
                                                                                                 coding_positions,
                                                                                                 coding_exon_starts)):
            counter += 1
            variant1.append(var1)
            variant2.append(var2)
            ppl = ppl_ids_w_pair(indf, var1, var2, cas)
            ppl_included_dict[counter] = inds.isin(ppl)
        else:
            continue

    targ_pairs_df = pd.DataFrame({'var1': variant1, 'var2': variant2})
    other_df = pd.DataFrame.from_dict(ppl_included_dict, orient='index')
    other_df.columns = inds
    out = pd.concat([targ_pairs_df, other_df], axis=1)

    return out


def main(args):
    gene = args['<gene>']
    filt = int(args['<filt>'])
    cas = args['<cas>']
    annots = ef.load_gene_annots(args['<annots>'])
    out_dir = args['<out_dir>']
    run_type = args['--type']
    geneid = ef.get_canonical(gene, annots, out_dir)
    indf = pd.read_hdf(args['<hap_targs>'], gene)
    if run_type == 'a':
        outdf = get_arcplot_df(indf, geneid, gene, cas, annots, out_dir).query('percent_pop_covered > @filt')
        outdf.to_csv(out_dir + f'shared_hets_{filt}.tsv', sep='\t', index=False)
    elif run_type == 's':
        outdf = get_set_cover_df(indf, geneid, gene, cas, annots, out_dir)
        outdf.to_csv(out_dir + f'set_cover_{gene}.tsv', sep='\t', index=False)


if __name__ == '__main__':
    arguments = docopt(__doc__)
    print(arguments)
    main(arguments)
