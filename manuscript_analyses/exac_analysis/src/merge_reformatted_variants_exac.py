import pandas as pd 
import os, sys

chrom=str(sys.argv[1])

store = pd.HDFStore(f'/pollard/data/projects/AlleleAnalyzer_data/exac_data/exac_formatted_variants_by_chrom/chr{chrom}_reformatted.h5')

region_files = pd.read_csv(f'/pollard/home/kathleen/projects/AlleleAnalyzer/manuscript_analyses/exac_analysis/dat/exac_regions_w_filenames.tsv',
	sep='\t', dtype={'chrom':str}).query('chrom == @chrom')
# region_files['chrom'] = region_files['chrom'].astype(str)
# region_files = region_files.query('chrom == @chrom')

for ix, row in region_files.iterrows():
	region = row['unique_name']
	f = row['gens_fname']
	region_df = pd.read_hdf(f)
	region_df['region'] = region
	min_item_sizes = {'chrom':2, 'pos':150, 'ref':500, 'alt':500, 'region':20}
	store.append('all', region_df, mode='a', index=False, format='table', append=True, 
		data_columns=['chrom','pos','region'], min_itemsize=min_item_sizes)

store.close()
print(f'Chromosome {chrom} complete.')