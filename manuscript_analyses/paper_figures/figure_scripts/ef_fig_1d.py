import pandas as pd
import sys
import numpy as np

"""
python ef_fig_1d.py single_targ_dfs/ dual_targ_dfs/ output_prefix
Kathleen Keough 2018
"""

# overall categories

dual_targetable = []
single_targetable = []
less_than_two_hets = []
no_hets = []
single_not_targ = []
dual_not_targ = []
no_coding_exons = []

# load single targ category lists

with open('wtc_hg19_single_targ/not_enough_hets.txt','r') as f:
    no_hets.extend(f.read().splitlines())
# with open('wtc_hg19_single_targ/no_targetable_inds.txt','r') as f:
    # single_not_targ.extend(f.read().splitlines())
with open('wtc_hg19_single_targ/no_coding_exons.txt','r') as f:
    no_coding_exons.extend(f.read().splitlines())
with open('wtc_hg19_single_targ/genes_evaluated.txt','r') as f:
    genes_eval_single = f.read().splitlines()

# load dual targ category lists

with open('wtc_hg19_whole_genome/not_enough_hets.txt','r') as f:
    less_than_two_hets.extend(f.read().splitlines())
with open('wtc_hg19_whole_genome/no_targetable_inds.txt','r') as f:
    dual_not_targ.extend(f.read().splitlines())
with open('wtc_hg19_whole_genome/no_coding_exons.txt','r') as f:
    no_coding_exons.extend(f.read().splitlines())
with open('wtc_hg19_whole_genome/genes_evaluated.txt','r') as f:
    genes_eval_dual = f.read().splitlines()

not_eval = []

# categorize each gene
for gene in genes_eval_single:
    try:
        df = pd.read_hdf(f'{sys.argv[1]}{gene}.h5')
        if df['targ_all'].any():
            single_targetable.append(gene)
        else:
            single_not_targ.append(gene)
    except FileNotFoundError:
        not_eval.append(gene)

for gene in genes_eval_dual:
    try:
        df = pd.read_hdf(f'{sys.argv[2]}{gene}.h5')
        if df['targ_all'].any():
            dual_targetable.append(gene)
        else:
            dual_not_targ.append(gene)
    except FileNotFoundError:
        not_eval.append(gene)

out_dict = {}
out_dict['Dual Targetable']=pd.DataFrame({'gene':list(set(dual_targetable))})
out_dict['Single Targetable']=pd.DataFrame({'gene':list(set(single_targetable))})
out_dict['Dual Not Targetable']=pd.DataFrame({'gene':list(set(dual_not_targ))})
out_dict['Single Not Targetable']=pd.DataFrame({'gene':list(set(single_not_targ))})
out_dict['No Coding Exons']=pd.DataFrame({'gene':list(set(no_coding_exons))})
out_dict['Less Than 2 Het Variants']=pd.DataFrame({'gene':list(set(less_than_two_hets))})
out_dict['No Het Variants']=pd.DataFrame({'gene':list(set(no_hets))})

out_prefix = sys.argv[3]
np.save(f'{out_prefix}categorized_variants.npy',out_dict)

with open(f'{out_prefix}not_eval.txt','w') as f:
    for gene in not_eval:
        f.write(gene + '\n')