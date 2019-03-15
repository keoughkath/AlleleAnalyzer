import pandas as pd
import sys

prefix = sys.argv[1]
out = sys.argv[2]

superpop_dict = {'ALL':'Total\n Population','EAS':'East\n Asian', 'SAS':'South\n Asian','AFR':'African','AMR':'Admixed\n American',
                'EUR':'European'}

dfs = []
for superpop in ['ALL','EAS','SAS','EUR','AMR','AFR']:
    indf = pd.read_csv(prefix + superpop + '.tsv', sep='\t')
    indf['pop'] = superpop_dict[superpop]
    dfs.append(indf)    

alldf = pd.concat(dfs)

alldf.to_csv(out + 'all_pops_arcplot_input.tsv', sep='\t', index=False)