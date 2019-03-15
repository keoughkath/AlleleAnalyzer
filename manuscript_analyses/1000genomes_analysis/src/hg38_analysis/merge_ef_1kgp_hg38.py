import pandas as pd 
import os, sys
from functools import reduce

chrom='chr' + str(sys.argv[1])

store = pd.HDFStore(f'/pollard/data/projects/AlleleAnalyzer_data/1kgp_data/hg38_analysis/1kgp_excisionfinder_results/results_by_chrom/{chrom}_ef_results.h5')

def translate_gene_name(gene_name):
    """
    HDF5 throws all sort of errors when you have weird punctuation in the gene name, so
    this translates it to a less offensive form.
    """
    repls = ("-", "dash"), (".", "period")
    trans_gene_name = reduce(lambda a, kv: a.replace(*kv), repls, str(gene_name))
    return trans_gene_name

with open('/pollard/home/kathleen/projects/AlleleAnalyzer/manuscript_analyses/1000genomes_analysis/get_gene_list/genes_hg38.txt','r') as f:
	genes = f.read().splitlines()

for gene in genes:
	fpath = f'/pollard/data/projects/AlleleAnalyzer_data/1kgp_data/hg38_analysis/1kgp_excisionfinder_results/results_by_chrom/{chrom}_ef_results/{gene}.h5'
	if os.path.exists(fpath):
		region_df = pd.read_hdf(fpath)
		gene_name_t = translate_gene_name(gene)
		store.append(gene_name_t, region_df, mode='a', index=False, format='table', append=True, 
			data_columns=['sample'])
	else:
		continue

store.close()
print(f'Chromosome {chrom} complete.')