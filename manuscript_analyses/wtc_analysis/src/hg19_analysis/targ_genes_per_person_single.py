"""
python targ_genes_per_person.py genes_evaluated.txt samples.txt targ_hdfs_dir/ out_prefix
Kathleen Keough 2018.
"""

import pandas as pd
import sys
import numpy as np

# load genes evaluated
with open(sys.argv[1],'r') as f:
    genes_eval = pd.read_csv(f, sep='\t')#f.read().splitlines()

# set up dict to catch targetable genes per person
targ_genes = {}

with open(sys.argv[2],'r') as f:
	people = f.read().splitlines()
	
for person in people:
    targ_genes[person] = []

not_eval = []

# get targetable genes per person
for ix, row in genes_eval.iterrows():
    gene = row['official_gene_symbol']
    chrom = row['chrom']
    try:
        df = pd.read_hdf(f'{sys.argv[3]}{chrom}_results/{gene}.h5')
        targ_ppl = df.query('targ_all')['sample'].tolist()
        for person in targ_ppl:
            targ_genes[person].append(gene)
    except FileNotFoundError:
        not_eval.append(gene)

out_prefix = sys.argv[4]

# save this dict to output
np.save(f'{out_prefix}genes_per_person.npy',targ_genes)

# save not evaluated
with open(f'{out_prefix}not_eval.txt','w') as f:
    for gene in not_eval:
        f.write(f'{gene}\n')