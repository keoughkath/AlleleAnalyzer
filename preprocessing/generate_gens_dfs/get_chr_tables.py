#!/usr/bin/env python3
# Rewriting `get_chr_tables.sh` in python
# -*- coding: utf-8 -*-
"""
get_chr_tables.py generates a table (tsv file) listing all variants in a defined interval for a specified 
individual (based on input VCF file). This basically reformats genotypes from VCF for easier 
processing later when designing sgRNAs.
Written in Python v 3.6.1.
Kathleen Keough and Michael Olvera 2017-2018.
Usage:
	get_chr_tables.py <vcf_file> <locus> <outdir> <name>

Arguments:
	vcf_file           The sample vcf file, separated by chromosome. 
	locus			   Locus from which to pull variants, in format chromosome:start-stop.
	outdir             Directory in which to save the output files.
	name			   The name for the output file.
"""
import pandas as pd
from docopt import docopt
import subprocess, os, sys
import regex as re

__version__='0.0.0'

REQUIRED_BCFTOOLS_VER = 1.5

def norm_chr(chrom_str):
	return chrom_str.replace('chr','')

def check_bcftools():
	""" 
	Checks bcftools version, and exits the program if the version is incorrect
	"""
	version = subprocess.run("bcftools -v | head -1 | cut -d ' ' -f2", shell=True,\
	 stdout=subprocess.PIPE).stdout.decode("utf-8").rstrip()
	if float(version) >= REQUIRED_BCFTOOLS_VER:
		print(f'bcftools version {version} running')

	else: 
		print(f"Error: bcftools must be >=1.5. Current version: {version}")
		exit()

def fix_multiallelics(cell):
	"""
	bcftools doesn't complete splitting of multiallelic variant sites.
	:param cell: genotype, str.
	:return: genotype as is if not multiallelic otherwise split multiallelic genotype, str.
	"""
	splitters = [',', ';']
	if any(splitter in str(cell) for splitter in splitters):
		cell = re.split(';|,', cell)[0]
	return cell


def het(genotype):
	gen1, gen2 = genotype.split('/|\|')
	return gen1 != gen2


def filter_hets(gens_df, mode='keep het'):
	"""
	if user specifies that they want homozygous instead, specify that here (not implemented yet)
	"""
	gens_df['het'] = gens_df.apply(lambda row: het(row['genotype']), axis=1)
	print(gens_df.head())
	return gens_df.query('het')[['chrom', 'pos', 'ref', 'alt', 'genotype']]


def main(args):
	# Make the outdir
	if not os.path.exists(args['<outdir>']):
		os.makedirs(args['<outdir>'])

	# Check if bcftools is installed, and then check version number
	check_bcftools()

	# See if chrom contains chr
	chrom = norm_chr(args['<locus>'].split(':')[0])
	start = args['<locus>'].split(':')[1].split('-')[0]
	stop = args['<locus>'].split(':')[1].split('-')[1]

	bcl_v=f"bcftools view -r chr{args['<chrom>']}:{str(start)}-{str(stop)} {args['<vcf_file>']}"
	
	# Pipe for bcftools
	bcl_view = subprocess.Popen(bcl_v,shell=True, stdout=subprocess.PIPE)
	bcl_norm = subprocess.Popen("bcftools norm -m -",shell=True, stdin=bcl_view.stdout, stdout=subprocess.PIPE)
	bcl_query = subprocess.Popen("bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n'",shell=True,
	 stdin=bcl_norm.stdout, stdout=subprocess.PIPE)
	bcl_query.wait() # Don't do anything else untill bcl_query is done running.

	# output  
	raw_dat = bcl_query.communicate()[0].decode("utf-8")

	temp_file_name=f"{args['<outdir>']}/{args['<chrom>']}_prechrtable.txt"
	with open(temp_file_name, 'w') as f:
		f.write(raw_dat)
		f.close()

	# Append fix_chr_tables.py
	vars = pd.read_csv(temp_file_name, sep='\t', header=None, names=['chrom', 'pos', 'ref', 'alt', 'genotype'],
		usecols=['chrom', 'pos', 'ref', 'alt', 'genotype'])

	if vars.empty:
		print('No heterozygous variants in this region for this individual. Exiting.')
		exit()

	if 'chr' in str(vars.chrom.iloc[0]):
		vars['chrom'] = vars['chrom'].map(lambda x: norm_chr(x))

	vars = vars.query('chrom == @chrom')
	vars_fixed = vars.applymap(fix_multiallelics)

	if args['<name>']:
		outname = f"{args['<name>']}.hdf5"
	else:
		outname = f'chr{chrom}_gens.hdf5'

	vars_fixed.to_hdf(os.path.join(args['<outdir>'], outname), 'all', complib='blosc')

	os.remove(temp_file_name)


if __name__ == '__main__':
	arguments = docopt(__doc__, version='0.1')
	main(arguments)
