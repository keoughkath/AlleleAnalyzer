# get gens for all chromosomes

for chrom in {1..22}; do
	echo $chrom
	bash ../excisionFinder/preprocessing/generate_gens_dfs/get_chr_tables_1kgp.sh /pollard/data/1kg/phase3/hg19/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf $chrom 1kgp_hg19/hg19_gens/
done