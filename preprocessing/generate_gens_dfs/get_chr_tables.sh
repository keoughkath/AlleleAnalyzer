# get chr tables
# required that you have bcftools installed, at least version 1.5

bcf_fn=$1
chrom=$2
outdir=$3

# The -p option will only write the directory if it doesn't exist
mkdir -p $outdir

echo $1 $2 $3 


# Check if bcftools is installed, and then check version number
version=`bcftools -v | head -1 | cut -d ' ' -f2`

req_version=1.5

if [ "$version" == "$req_version" ]; then
	echo "bcftools version $version, running"
else
	echo "Error: bcftools must be 1.5. Current version: $version"
	exit 1
fi

sample_list=`bcftools query -l $bcf_fn`

bcftools norm -m - ${bcf_fn} | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' --samples $sample_list > $outdir/${chrom}_prechrtable.txt

python fix_chr_tables.py $outdir/${chrom}_prechrtable.txt $chrom $outdir

rm $outdir/${chrom}_prechrtable.txt
