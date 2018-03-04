import pandas as pd
import sys
import crisprtools as cpt

infile = sys.argv[1]
gene_vars_file = sys.argv[2]
outfile = sys.argv[3]
track_name = sys.argv[4].replace(' ','_')

gene_bed = pd.read_csv(infile, sep='\t')
gene_vars = pd.read_csv(gene_vars_file, sep='\t',header=None,
                       names=['chrom','variant_position','ref','alt','rsID','AF'])

def adjusted_length(row):
	"""
	Adds on the length of the PAM to the sequnce length.
	"""
	cas = row['cas_type']
	if row['strand'] == '+':
		return (row['start'],row['stop']+cpt.tpp_for[cas][1]) 
	else: 
		return (row['start']-cpt.tpp_for[cas][1],row['stop'])


gene_bed['full_start'], gene_bed['full_stop'] = zip(* gene_bed.apply(adjusted_length, axis=1))

gene_vars['chrom'] = list(map(lambda x: 'chr' + str(x), gene_vars['chrom']))
gene_bed['variant position in guide'] = gene_bed['variant_position_in_guide']
#gene_bed['score'] = 1000*(1/(gene_bed['variant_position_in_guide'] + 1))
gene_bed = gene_bed.merge(gene_vars, how='left')
gene_bed['label'] = gene_bed.apply(lambda row: row['guide_id']+'_'+str(row['variant_position_in_guide']), axis=1) 

gene_bed_display = gene_bed.query('variant_position_in_guide != 2 and variant_position_in_guide >= 0')[['chrom',
                                                                     'full_start','full_stop',
                                                                     'label','norm_score',
                                                                    'strand','start','stop']]

header_str = f'track name={track_name} description=AS cut sites as produced by ExcisionFinder visibility=3 useScore=1 '

gene_bed_display.to_csv(outfile,sep='\t',index=False,
                       header=[header_str,'','','','','','',''])