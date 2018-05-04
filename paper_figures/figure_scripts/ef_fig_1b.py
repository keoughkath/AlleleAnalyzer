import pandas as pd

cas_list=['SpCas9','SpCas9_VRER','SpCas9_EQR','SpCas9_VQR_1','SpCas9_VQR_2','StCas9','StCas9_2','SaCas9','SaCas9_KKH','nmCas9','cjCas9','cpf1']

in_dict = {}
near_dict = {}
both_dict = {}

for cas in cas_list:
	in_dict[cas] = 0
	near_dict[cas] = 0
	both_dict[cas] = 0

global total_vars
total_vars = 0

for chrom in list(range(1,23)):
	chrom=str(chrom)
	print(chrom)
	annot_variants=pd.read_hdf(f'/pollard/home/kathleen/projects/genome_surgery/1kgp_hg19/hg19_targs/{chrom}_targ.hdf5','all').drop(columns=['chrom','pos','ref','alt'])
	annot_variants['id'] = annot_variants.index
	# get # of variants in chromosome annotated
	n_vars = annot_variants.shape[0]
	total_vars += n_vars
	for cas in cas_list:
		if cas == 'SpCas9_VQR_1':
			in_pam=set(annot_variants.query(f'makes_{cas} or breaks_{cas} or makes_SpCas9_VQR_2 or breaks_SpCas9_VQR_2')['id'].tolist())
			near_pam=set(annot_variants.query(f'var_near_{cas} or var_near_SpCas9_VQR_2')['id'].tolist())
			both = in_pam.intersection(near_pam)
			both_dict['SpCas9_VQR'] += len(both)
			in_only = in_pam.difference(near_pam)
			in_dict['SpCas9_VQR'] += len(in_only)
			near_only = near_pam.difference(in_pam)
			near_dict['SpCas9_VQR'] += len(near_only)
		elif cas == 'SpCas9_VQR_2':
			continue 
		else:
			in_pam=set(annot_variants.query(f'makes_{cas} or breaks_{cas}')['id'].tolist())
			near_pam=set(annot_variants.query(f'var_near_{cas}')['id'].tolist())
			both = in_pam.intersection(near_pam)
			both_dict[cas] += len(both)
			in_only = in_pam.difference(near_pam)
			in_dict[cas] += len(in_only)
			near_only = near_pam.difference(in_pam)
			near_dict[cas] += len(near_only)

in_df = pd.DataFrame.from_dict(in_dict, orient='index')
in_df.columns = ['in_pam']

near_df = pd.DataFrame.from_dict(near_dict, orient='index')
near_df.columns = ['near_pam']

both_df = pd.DataFrame.from_dict(both_dict, orient='index')
both_df.columns = ['both']

# set up output dataframe

plot_df_out = in_df.merge(near_df,
	left_index=True, right_index=True).merge(both_df, left_index=True, 
	right_index=True).divide(total_vars)

# save to file
plot_df_out.to_csv('figure_data/vars_near_in_df.tsv', sep='\t')




