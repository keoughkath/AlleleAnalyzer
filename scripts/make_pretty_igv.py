import pandas as pd
import sys

infile = sys.argv[0]
gene_vars_file = sys.argv[1]
outfile = sys.argv[2]

gene_bed = pd.read_csv(infile, sep='\t')
gene_vars = pd.read_csv(gene_vars_file, sep='\t',header=None,
                       names=['chrom','variant_position','ref','alt','rsID','AF'])

gene_vars['chrom'] = list(map(lambda x: 'chr' + str(x), gene_vars['chrom']))
gene_bed['variant position in guide'] = gene_bed['variant_position_in_guide']
gene_bed['score'] = 1000*(1/(gene_bed['variant_position_in_guide'] + 1))
gene_bed = gene_bed.merge(gene_vars, how='left')

gene_bed_display = gene_bed.query('variant_position_in_guide != 2 and variant_position_in_guide >= 0')[['chrom',
                                                                     'start','stop',
                                                                     'variant position in guide','score',
                                                                    'strand']]

gene_bed_display.to_csv(outfile,sep='\t',index=False,
                       header=['track name=gene_AS description=AS cut sites gene WTC visibility=3 useScore=1 ','','','','',''])