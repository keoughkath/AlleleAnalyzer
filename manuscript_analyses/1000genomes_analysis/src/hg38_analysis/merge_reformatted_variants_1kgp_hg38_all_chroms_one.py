import pandas as pd 
import os

store = pd.HDFStore(f'/pollard/data/projects/AlleleAnalyzer_data/1kgp_data/hg38_analysis/1kgp_formatted_variants_by_chrom/1kgp_hg38_formatted_variants_all.h5',
	complib='blosc', complevel=9)

for chrom in list(range(1,23)) + ['X','Y']:
	df = pd.read_hdf(f'/pollard/data/projects/AlleleAnalyzer_data/1kgp_data/hg38_analysis/1kgp_formatted_variants_by_chrom/chr{chrom}_reformatted.h5')
	chrom = df.head()['chrom'].tolist()[0]
	min_item_sizes = {'chrom':4, 'pos':150, 'ref':500, 'alt':500, 'region':20}
	store.append('chr' + str(chrom), df, mode='a', index=False, format='table', append=True, 
		data_columns=['chrom','pos','region'], min_itemsize=min_item_sizes)
	print(f'Chromosome {chrom} complete.')

store.close()
print(f'All chromosomes complete.')