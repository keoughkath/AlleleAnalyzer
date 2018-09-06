# get PAM positions for each chromosome in hg19

import crisprtools
import numpy as np
import sys
from pyfaidx import Fasta

# keep track of chromosome since this will be run with bash script

chrom = sys.argv[1]

# get access to hg19

hg19 = Fasta('/pollard/home/kathleen/hg19/hg19.fa',as_raw=True)

cas_list = ['SpCas9','SpCas9_VRER','SpCas9_EQR','SpCas9_VQR_1',
'SpCas9_VQR_2','StCas9','StCas9_2','SaCas9','SaCas9_KKH','nmCas9','cjCas9']

# get set of positions for each type of cas

for cas in cas_list:
	for_starts, rev_starts = crisprtools.find_spec_pams(cas,str(hg19[str(chrom)]))
	savestr_for = 'hg19_chr_pams_june_2018/'+str(chrom)+'_'+str(cas) + '_pam_sites_for.npy'
	savestr_rev = 'hg19_chr_pams_june_2018/'+str(chrom)+'_'+str(cas) + '_pam_sites_rev.npy'
	np.save(savestr_for,list(for_starts))
	np.save(savestr_rev,list(rev_starts))

# cpf1 is 5pp

for_starts, rev_starts = crisprtools.find_spec_pams('cpf1',str(hg19[str(chrom)]),orient='5prime')
np.save('hg19_chr_pams_june_2018/'+str(chrom)+'_'+'cpf1_pam_sites_for.npy',list(for_starts))
np.save('hg19_chr_pams_june_2018/'+str(chrom)+'_'+'cpf1_pam_sites_rev.npy',list(rev_starts))