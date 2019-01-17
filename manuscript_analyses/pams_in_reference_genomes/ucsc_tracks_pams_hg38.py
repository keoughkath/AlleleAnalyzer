# make UCSC genome browser track

import pandas as pd

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

pams = pd.read_csv('/pollard/data/projects/AlleleAnalyzer_data/pam_beds_hg38/pams_hg38.bed',
                   sep='\t', header=None, skiprows=1, names=['chrom','start','stop','cas','blah','strand',
                                                            'start_','stop_','color'])
pams['cas_'] = pams['cas'].str[:-1]
pams['color_'] = pams['cas_'].map(cas_to_colors)

for cas in cas_to_colors.keys():
    print(cas)
    cas_color = cas_to_colors[cas]
    pams.query('cas_ == @cas')[['chrom','start','stop','cas_','blah','strand',
      'start_','stop_','color_']].to_csv(f'/pollard/data/projects/AlleleAnalyzer_data/pam_beds_hg38/ucsc_browser/{cas}_ucsc.bed',
           sep='\t', index=None, header=None)