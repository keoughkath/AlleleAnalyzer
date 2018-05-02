import pandas as pd

for chrom in list(range(1,23)):
	chrom = str(chrom)
	annot_variants = pd.read_hdf(f')