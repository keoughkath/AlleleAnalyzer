import pandas as pd 
import os, sys

chrom='chr' + str(sys.argv[1])

store = pd.HDFStore(f'/pollard/data/projects/AlleleAnalyzer_data/1kgp_data/hg19_analysis/1kgp_formatted_variants_by_chrom/{chrom}_reformatted.h5')

region_files = pd.read_csv(f'/pollard/data/projects/AlleleAnalyzer_data/1kgp_data/hg19_analysis/region_annots.tsv',
	sep='\t').query('chrom == @chrom')

for ix, row in region_files.iterrows():
	region = row['region_id']
	f = row['gens_fname']
	region_df = pd.read_hdf(f)
	region_df['region'] = region
	min_item_sizes = {'chrom':2, 'pos':150, 'ref':500, 'alt':500, 'region':20}
	store.append('all', region_df, mode='a', index=False, format='table', append=True, 
		data_columns=['chrom','pos','region'], min_itemsize=min_item_sizes)

store.close()
print(f'Chromosome {chrom} complete.')