"""
Get number of variants annotated from 1KGP analysis.
Kathleen Keough 2018.

python n_vars_annotated targ_dfs_dir/ out_prefix
"""

import pandas as pd
from io import StringIO
import multiprocessing as mp
import glob
import sys

n_vars = {}

def get_n_vars(chrom):
    df = pd.read_hdf(f'{sys.argv[1]}/{chrom}_targ.hdf5','all')
    n_vars = df.shape[0]
    n_vars[chrom] = n_vars

count = min(mp.cpu_count(),22) # only use as many CPUs as available
pool = mp.Pool(processes=count)
chroms = list(range(1,23,1))

pool.map(get_n_vars, chroms)
pool.close()
pool.join()

n_vars_df = pd.DataFrame.from_dict(n_vars, orient='index')
n_vars_df.to_csv(f'{sys.argv[2]}_nvars.tsv', sep='\t')