import pandas as pd 
import os

store = pd.HDFStore(f'/pollard/data/projects/AlleleAnalyzer_data/exac_data/exac_annotated_variants_by_chrom/exac_annots_all.h5')

for chrom in list(range(1,23)) + ['X','Y']:
	df = pd.read_hdf(f'/pollard/data/projects/AlleleAnalyzer_data/exac_data/exac_annotated_variants_by_chrom/chr{chrom}_annots.h5')
	chrom = str(df.head()['chrom'].tolist()[0])
	min_item_sizes = {**{'chrom':2, 'pos':150, 'region':150, 'ref':500, 
	'alt':500}, **dict(zip(df.columns[5:], [250] * len(df.columns[5:])))}
	store.append('chr' + str(chrom), df, mode='a', index=False, format='table', append=True, 
		data_columns=['chrom','pos','region'], min_itemsize=min_item_sizes)
	print(f'Chromosome {chrom} complete.')

store.close()
print(f'All chromosomes complete.')