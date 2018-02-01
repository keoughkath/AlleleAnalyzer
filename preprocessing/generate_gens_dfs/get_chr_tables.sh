# get chr tables
# required that you have bcftools installed, at least version 1.5

bcf_fn=$1
chrom=$2
outdir=$3

mkdir $outdir

echo $1 $2 $3 

sample_list=`~/tools/bcftools-1.5/./bcftools query -l $bcf_fn`

~/tools/bcftools-1.5/./bcftools norm -m - ${bcf_fn} | ~/tools/bcftools-1.5/./bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' --samples $sample_list > $outdir/${chrom}_prechrtable.txt

python fix_chr_tables.py $outdir/${chrom}_prechrtable.txt $chrom $outdir

rm $outdir/${chrom}_prechrtable.txt
