# get gens for all chromosomes

for chrom in {1..22}; do
	echo $chrom
	bash ../excisionFinder/preprocessing/generate_gens_dfs/get_chr_tables_1kgp.sh /pollard/data/1kg/phase3/hg38/ALL.chr${chrom}_GRCh38.genotypes.20170504.bcf $chrom 1kgp_hg38/hg38_gens/
done