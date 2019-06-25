# !/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script takes in the long-form data that ExcisionFinder puts out in "exhaustive" mode and converts it into a format that can be used by arcplot_generic.R.

Usage: 
	gen_arcplot_input.py [-h] <input_df> <out> [--sample_legend=<sample_legend>] [--pop=<population>]

Arguments:
	input_df   			Long-form input df outputted by ExcisionFinder in "exhaustive" mode.
	out   				Prefix for output file.

Options:
	-h   			    Shows this page and exit.
	--pop=<population>  Analyze a specific population, e.g. AMR.
	--sample_legend=<sample_legend>       File with population data from 1KGP. Located here: http://lighthouse.ucsf.edu/public_files_no_password/excisionFinderData_public/1kgp_dat/
"""

import pandas as pd
from docopt import docopt

__version__ = '0.0.0'

def filt_pops(df, sample_legend, pop='all'):
    df = df.merge(sample_legend, how='left')
    if pop==None:
        return df, sample_legend.shape[0]
    else:
        return df.query('pop == @pop or superpop == @pop'), sample_legend.query('superpop == @pop or pop == @pop').shape[0]

def main(args):
	df = pd.read_csv(args['<input_df>'], sep='\t', low_memory=False)
	if args['--sample_legend']:
		sample_legend = pd.read_csv(args['--sample_legend'], sep='\t', header=0,
		                           names=['superpop','pop','sex'])
		sample_legend['ind'] = sample_legend.index

		all_pops, n_pop = filt_pops(df, sample_legend, args['--pop'])
	else:
		all_pops = df
		n_pop = int(input('Enter number of individuals: '))

	inds_per_pair = all_pops[['var1','var2','ind']].drop_duplicates().groupby(['var1','var2']).ind.count().reset_index()
	inds_per_pair['perc_shared'] = inds_per_pair['ind'] / n_pop * 100.0
	inds_per_pair.columns = ['var1','var2','n_inds','percent_pop_covered']

	inds_per_pair.to_csv(args['<out>'] + '.tsv', sep='\t', index=False)

if __name__ == '__main__':
    arguments = docopt(__doc__, version=__version__)
    print(arguments)
    main(arguments)