import pandas as pd 
import os

store = pd.HDFStore(f'/pollard/data/projects/AlleleAnalyzer_data/1kgp_data/hg19_analysis/1kgp_annotated_variants_by_chrom/1kgp_hg19_annots_all.h5',
	complib='blosc', complevel=9)

for chrom in list(range(1,23)) + ['X','Y']:
	df = pd.read_hdf(f'/pollard/data/projects/AlleleAnalyzer_data/1kgp_data/hg19_analysis/1kgp_annotated_variants_by_chrom/chr{chrom}_annotated.h5')
	chrom = str(df.head()['chrom'].tolist()[0])
	min_item_sizes = {**{'chrom':2, 'pos':150, 'region':150, 'ref':500, 
	'alt':500}, **dict(zip(df.columns[5:], [250] * len(df.columns[5:])))}
	store.append('chr' + str(chrom), df, mode='a', index=False, format='table', append=True, 
		data_columns=['chrom','pos','region'], min_itemsize=min_item_sizes)
	print(f'Chromosome {chrom} complete.')

store.close()
print(f'All chromosomes complete.')