for chrom in {1..22} X Y
do
	echo ${chrom}
bcftools view -H /pollard/home/kathleen/conklin_wt_seq_data/wtc_wgs_data/phased_yin/wtc_PASS_hg19.phased.vcf.gz | cut -f1 | head -1 | bcftools view -g ^miss -r chr${chrom} /pollard/home/kathleen/conklin_wt_seq_data/wtc_wgs_data/phased_yin/wtc_PASS_hg19.phased.vcf.gz | bcftools norm -m - | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\n' > /pollard/data/projects/AlleleAnalyzer_data/wtc_data/hg19/wtc_formatted_variants_by_chrom2/chr${chrom}.tsv
done