import pandas as pd
import numpy as np
import sys
from functools import reduce

"""
ef_fig_1e.py genes_eval_file targ_h5_location out
Kathleen Keough 2018
"""

cas_list=['SpCas9','SpCas9_VRER','SpCas9_EQR','SpCas9_VQR_1','SpCas9_VQR_2','StCas9','StCas9_2','SaCas9','SaCas9_KKH','nmCas9','cjCas9','cpf1']

# catch output

cases = []
genes = []
perc_ppl_targ = []

genes_eval_file = sys.argv[1]

# get evaluated genes list
# with open(genes_eval_file,'r') as f:
#     genes_eval = f.read().splitlines()

def translate_gene_name(gene_name):
    """
    HDF5 throws all sort of errors when you have weird punctuation in the gene name, so
    this translates it to a less offensive form.
    """
    repls = ("-", "dash"), (".", "period")
    trans_gene_name = reduce(lambda a, kv: a.replace(*kv), repls, str(gene_name))
    return trans_gene_name

genes_eval = pd.read_csv(genes_eval_file, sep='\t')

n_genes = len(genes_eval)
counter = -1
not_eval = []
for ix, row in genes_eval.iterrows():
    gene = row['official_gene_symbol']
    chrom = row['chrom']
    print(gene, chrom)
    counter += 1
    translated_gene_name = translate_gene_name(gene)
    try:
        print(n_genes - counter)
        df = pd.read_hdf(f'{sys.argv[2]}/{chrom}_ef_results/{gene}.h5')
        genes.append(gene)
        cases.append('all')
        perc_ppl_targ.append(df.query('targ_all').shape[0]/2504.0)
        for cas in cas_list:
            if cas == 'SpCas9_VQR_1':
                genes.append(gene)
                cases.append('SpCas9_VQR')
                perc_ppl_targ.append(df.query('targ_SpCas9_VQR_1 or targ_SpCas9_VQR_2').shape[0]/2504.0)
            elif cas == 'SpCas9_VQR_2':
                continue
            else:
                genes.append(gene)
                cases.append(cas)
                perc_ppl_targ.append(df.query(f'targ_{cas}').shape[0]/2504.0)
    except FileNotFoundError:
        not_eval.append(gene)
        pass
    
# aggregate into df

overall_df = pd.DataFrame({'Cas':cases,'Gene':genes,
    '% people targetable':perc_ppl_targ})

with open(sys.argv[3] + 'not_eval.txt','w') as f:
    for gene in not_eval:
        f.write(gene+'\n')
overall_df.to_csv(sys.argv[3] + 'targ_per_gene_and_cas.tsv', sep='\t', index=False)