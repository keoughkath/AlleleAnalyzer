# bash script to find PAM sites in all chromosomes in a fasta

export FASTA=$1
export NCHROMS=$2
export OUTDIR=$3

for CHROM in $(seq 1 $NCHROMS)
do
	echo $CHROM
	python pam_pos_genome.py chr$CHROM $FASTA cpf1,SpCas9,SpCas9_VRER,SpCas9_EQR,SpCas9_VQR_1,SpCas9_VQR_2,StCas9,StCas9_2,SaCas9,SaCas9_KKH,nmCas9,cjCas9 $OUTDIR/
done