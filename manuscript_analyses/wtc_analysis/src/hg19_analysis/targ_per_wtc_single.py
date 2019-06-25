import pandas as pd
import numpy as np

targ = []
not_targ = []
not_eval = []

genes = pd.read_csv('/pollard/home/kathleen/projects/AlleleAnalyzer/manuscript_analyses/1000genomes_analysis/get_gene_list/gene_list_hg19.tsv',
                   sep='\t')

for ix, row in genes.iterrows():
    chrom = row['chrom']
    gene = row['official_gene_symbol']
    try:
        df = pd.read_hdf(f'/pollard/data/projects/AlleleAnalyzer_data/wtc_data/hg19/wtc_single_targ2/{chrom}_results/{gene}.h5')
        if df['targ_all'].item():
            targ.append(gene)
        else:
            not_targ.append(gene)
    except FileNotFoundError:
        not_eval.append(gene)

targ_genes = {}
targ_genes['WTC'] = targ


np.save('/pollard/data/projects/AlleleAnalyzer_data/wtc_data/hg19/wtc_targ_genes_per_person/single_targgenes_per_person.npy', targ_genes)

