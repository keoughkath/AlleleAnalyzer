# get gens for all chromosomes

for chrom in {1..22}; do
	echo $chrom
	bash get_chr_tables_wtc.sh wtc_phased_hg19.bcf $chrom wtc_hg19_whole_genome/hg19_gens/
done