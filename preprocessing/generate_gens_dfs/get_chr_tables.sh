# get chr tables
# required that you have bcftools installed, at least version 1.5

bcf_fn=$1
chrom="${2//chr}" # removes the 'chr', if included
outdir=$3

# The -p option will only write the directory if it doesn't exist
mkdir -p $outdir

echo $1 $2 $3 


# Check if bcftools is installed, and then check version number
version=`bcftools -v | head -1 | cut -d ' ' -f2`

req_version=1.5
# Checks to see if the version is greater than or equal to 1.5.
if [ "$(printf '%s\n' "$req_version" "$version" | sort | head -n1)" = "$req_version" ]; then
	echo "bcftools version $version, running"
else
	echo "Error: bcftools must be >=1.5. Current version: $version"
	exit 1
fi

sample_list=`bcftools query -l $bcf_fn`

# Removed --samples in the 'bcftools query' command, it would not let me run with multivariant vcfs.
bcftools norm -m - ${bcf_fn} | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n'  > $outdir/${chrom}_prechrtable.txt

# Added $(dirname "$0"), so that is the script is called from another dir, this script can be accessed
# Let me know what you think about this. Personally I like to run scripts from the directory where all my files are, or even
# from a higher directory that contains both my scripts and my input/output files. I got the sample output whether
# I ran the scrip from 'generate_gens_dfs' or not.
python $(dirname "$0")/fix_chr_tables.py $outdir/${chrom}_prechrtable.txt $chrom $outdir

rm $outdir/${chrom}_prechrtable.txt
