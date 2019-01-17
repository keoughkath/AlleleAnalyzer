import pandas as pd
import numpy as np
import gzip
import sys
import os
from pyfaidx import Fasta
import matplotlib as mpl
import matplotlib.pyplot as plt

# Get absolute path for gen_targ_dfs.py, and edit it for cas_object.py
ef_path = os.path.dirname(os.getcwd())
cas_obj_path = ef_path.replace('manuscript_analyses','scripts/')
sys.path.append(cas_obj_path)
import cas_object as cas_obj

# load PAM sites for one chromosome to start

cas_list=['SpCas9','SpCas9_VRER','SpCas9_EQR','SpCas9_VQR_1',
          'SpCas9_VQR_2','StCas9','StCas9_2','SaCas9',
          'SaCas9_KKH','nmCas9','cjCas9','cpf1']

chroms = list(range(1,23)) + ['Y'] + ['X']

def get_seq(row):
    seq = hg19['chr' + str(chrom)][row['pam_start']:row['pam_stop']]
    return seq

# map Cas proteins to colors

cas_to_colors = {
'SpCas9':'0,0,255',
'SpCas9_VRER':'102,0,204',
'SpCas9_EQR':'255,0,255',
'SpCas9_VQR_1':'255,0,102',
'SpCas9_VQR_2':'255,0,0',
'StCas9':'255,102,0',
'StCas9_2':'255,204,0',
'SaCas9':'0,153,0',
'SaCas9_KKH':'0,255,0',
'nmCas9':'0,204,153',
'cjCas9':'0,102,102',
'cpf1':'0,153,255'
}

bed_dfs = []
for chrom in chroms:
    for cas in cas_list:
        print(cas)
        cas_info = cas_obj.get_cas_enzyme(cas, os.path.join(cas_obj_path,'CAS_LIST.txt'))
        pam_size = len(cas_info.forwardPam)
        for_pams = []
        rev_pams = []
        chroms_for = []
        chroms_rev = []
        with gzip.open(f'/pollard/data/projects/AlleleAnalyzer_data/pam_sites_hg19/pam_sites_hg19_txt/chr{chrom}_{cas}_pam_sites_for.txt.gz', 'rb') as f:
            for line in f.readlines():
                for_pams.append(int(float(line.strip())) - 1)
                chroms_for.append('chr' + str(chrom))
        with gzip.open(f'/pollard/data/projects/AlleleAnalyzer_data/pam_sites_hg19/pam_sites_hg19_txt/chr{chrom}_{cas}_pam_sites_rev.txt.gz', 'rb') as f:
            for line in f.readlines():
                rev_pams.append(int(float(line.strip())))
                chroms_rev.append('chr' + str(chrom))
        if cas_info.primeness == "3'":
            bed_df_for = pd.DataFrame({'chrom':chroms_for,
                                      'pam_start':for_pams})
            bed_df_for['pam_stop'] = bed_df_for['pam_start'] + pam_size
            bed_df_rev = pd.DataFrame({'chrom':chroms_rev,
                                      'pam_stop':rev_pams})
            bed_df_rev['pam_start'] = bed_df_rev['pam_stop'] - pam_size
        else:
            bed_df_for = pd.DataFrame({'chrom':chroms_for,
                                      'pam_stop':for_pams})
            bed_df_for['pam_stop'] = bed_df_for['pam_stop'] - 1
            bed_df_for['pam_start'] = bed_df_for['pam_stop'] - pam_size
            bed_df_rev = pd.DataFrame({'chrom':chroms_rev,
                                      'pam_start':rev_pams})
            bed_df_rev['pam_start'] = bed_df_rev['pam_start'] - 1
            bed_df_rev['pam_stop'] = bed_df_rev['pam_start'] + pam_size
        bed_df_for['strand'] = '+'
        bed_df_rev['strand'] = '-'
        overall = pd.concat([bed_df_for, bed_df_rev], sort=False)
        overall['Cas'] = cas
        bed_dfs.append(overall)

bed_df = pd.concat(bed_dfs, sort=False)

bed_df['name'] = bed_df['Cas'] + bed_df['strand']
bed_df['score'] = 0
bed_df['thickStart'] = bed_df['pam_start']
bed_df['thickEnd'] = bed_df['pam_stop']
bed_df['itemRgb'] = bed_df['Cas'].replace(cas_to_colors)

bed_df[['chrom','pam_start','pam_stop','name','score','strand','thickStart',
               'thickEnd','itemRgb']].to_csv(f'/pollard/data/projects/AlleleAnalyzer_data/pam_beds_hg19/pams_hg19.bed', sep='\t', header=False, index=False)