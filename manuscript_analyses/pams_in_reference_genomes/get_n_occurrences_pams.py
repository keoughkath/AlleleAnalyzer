"""
Sum up # of occurrences of PAM site in genome for multiple varieties of Cas.
python get_n_occurrences_pams.py pams_dir/ out_prefix
Kathleen Keough 2018.
"""

import pandas as pd
import numpy as np
import sys

cas_list=['SpCas9','SpCas9_VRER','SpCas9_EQR','SpCas9_VQR_1','SpCas9_VQR_2','StCas9','StCas9_2','SaCas9','SaCas9_KKH','nmCas9','cjCas9','cpf1']
chrom_list = list(range(1,23)) + ['Y'] + ['X']

n_pams_dict = {}
n_pams_dict['SpCas9_VQR'] = 0

for cas in cas_list:
	print(cas)
	n_pams_dict[cas] = 0
	for chrom in chrom_list:
		print(chrom)
		for_pams = np.load(f'{sys.argv[1]}chr{chrom}_{cas}_pam_sites_for.npy')
		rev_pams = np.load(f'{sys.argv[1]}chr{chrom}_{cas}_pam_sites_rev.npy')
		if cas == 'SpCas9_VQR_1' or cas == 'SpCas9_VQR_2':
			n_pams_dict['SpCas9_VQR'] += len(for_pams) + len(rev_pams)
		else:
			n_pams_dict[cas] += len(list(for_pams)) + len(list(rev_pams))

np.save(f'{sys.argv[2]}n_pams_per_cas.npy',n_pams_dict)


