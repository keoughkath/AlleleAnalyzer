import pandas as pd 
import os, sys

chrom=sys.argv[1]

store = pd.HDFStore(f'/pollard/data/projects/AlleleAnalyzer_data/exac_data/exac_annotated_variants_by_chrom_parallel/chr{chrom}_annots.h5')

region_files = pd.read_csv(f'/pollard/data/projects/AlleleAnalyzer_data/exac_data/exac_regions_by_chrom/chr{chrom}_regions.tsv',
	sep='\t')

for ix, row in region_files.iterrows():
	region = row['unique_name']
	region_df = pd.read_hdf(row['annots_fname'])
	region_df['region'] = region
	min_item_sizes = {**{'chrom':2, 'pos':150, 'region':150, 'ref':500, 
	'alt':500}, **dict(zip(region_df.columns[5:], [250] * len(region_df.columns[5:])))}
	store.append('all', region_df, mode='a', index=False, format='table', append=True, 
		data_columns=['chrom','pos','region'], min_itemsize=min_item_sizes)

store.close()
print(f'Chromosome {chrom} complete.')