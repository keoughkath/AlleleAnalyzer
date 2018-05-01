# bash script to find all PAM sites (except cpf1) in hg19

for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
	echo $chrom
	python pam_pos_hg38.py $chrom
done