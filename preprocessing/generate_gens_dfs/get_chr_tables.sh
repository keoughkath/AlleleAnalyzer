# get chr tables
# required that you have bcftools installed, at least version 1.5

bcf_fn=$1
chrom="${2//chr}" # removes the 'chr', if included
outdir=$3
start=$4
stop=$5
name=$6

echo $bcf_fn
# The -p option will only write the directory if it doesn't exist
mkdir -p $outdir


# Check if bcftools is installed, and then check version number
version=`bcftools -v | head -1 | cut -d ' ' -f2`

# Checks to see if the version is greater than or equal to the required version.
req_version=1.5
if [ "$(printf '%s\n' "$req_version" "$version" | sort | head -n1)" = "$req_version" ]; then
	echo "bcftools version $version, running"
else
	echo "Error: bcftools must be >=1.5. Current version: $version"
	exit 1
fi

bcftools view -r chr${chrom}:${start}-${stop} ${bcf_fn} | bcftools norm -m - | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' > $outdir/${chrom}_prechrtable.txt

# Added $(dirname "$0"), so that is the script is called from another dir, this script can be accessed
python $(dirname "$0")/fix_chr_tables.py $outdir/${chrom}_prechrtable.txt $chrom $outdir $name

rm $outdir/${chrom}_prechrtable.txt