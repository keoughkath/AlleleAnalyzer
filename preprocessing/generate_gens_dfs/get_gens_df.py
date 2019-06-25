#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
get_gens_dfs.py generates a table (tsv file) listing all variants in a defined interval for a specified 
individual (based on input VCF file). This basically reformats genotypes from VCF for easier 
processing later when designing sgRNAs.
Written in Python v 3.6.1.
Kathleen Keough et al 2017-2018.

Usage:
    get_gens_dfs.py <vcf_file> <locus> <out> [-v] [--bed] [--chrom]

Arguments:
    vcf_file           The sample vcf file, separated by chromosome. BCF also supported. 
    locus               Locus from which to pull variants, in format chromosome:start-stop, or a BED file if --bed.
    out                   The name for the output file and directory in which to save the output files.
Options:
    -v                    Verbose mode.
    --bed              Indicates that a BED file is being used in place of a locus.
    --chrom            Run on entire chromosome.
"""

import pandas as pd
from docopt import docopt
import subprocess, os, sys, logging
import regex as re
from io import StringIO

# Append path to metadata script
ef_path = os.path.dirname(os.path.realpath(__file__))
metadata_path = ef_path.replace('generate_gens_dfs','')

sys.path.append(metadata_path)

from get_metadata import add_metadata

__version__='1.0.0'

REQUIRED_BCFTOOLS_VER = 1.5

def norm_chr(chrom_str, vcf_chrom):
    chrom_str = str(chrom_str)
    if not vcf_chrom:
        return chrom_str.replace('chr','')
    elif vcf_chrom and not chrom_str.startswith('chr'):
        return 'chr' + chrom_str
    else:
        return chrom_str

def check_bcftools():
    """ 
    Checks bcftools version, and exits the program if the version is incorrect
    """
    try:
        version = subprocess.run("bcftools -v | head -1 | cut -d ' ' -f2", shell=True,\
         stdout=subprocess.PIPE).stdout.decode("utf-8").rstrip()
        print(version)
        version = float(version.strip())
        # if float(version) >= REQUIRED_BCFTOOLS_VER:
        #     logging.info(f'bcftools version {version} running')
        # else: 
        #     logging.error(f"Error: bcftools must be >=1.5. Current version: {version}")
        #     exit()
    except ValueError:
        print('Error: Check that bcftools version >= 1.5 is installed. Exiting.')
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

def fix_natural_language(name):
    """
    Fixes NaturalNameWarning given by trying to write an hdf5 column name
    """
    for ch in r"\`*{}[]()>#+-.!$":
        if ch in name:
            name = name.replace(ch,"_")
    return name

def main(args):

    logging.info(args)
    vcf_in = args['<vcf_file>']
    out = args['<out>']
    # Check if bcftools is installed, and then check version number
    check_bcftools()

    # analyze regions specified in BED file
    if (args['<locus>'].endswith('.bed') and not args['--bed']) or (args['<locus>'].endswith('.BED') and not args['--bed']):
        logging.error('Must specify --bed if inputting a BED file. Exiting.')
        exit()
    elif args['--bed']:
        bed_file = args['<locus>']
        logging.info(f'Analyzing BED file {bed_file}')
        bed_df = pd.read_csv(bed_file, sep='\t', header=None, comment='#', names=['chr','start','stop','locus'])
        vcf_chrom = subprocess.Popen(f'bcftools view -H {vcf_in} | cut -f1 | head -1', shell=True, 
            stdout=subprocess.PIPE).communicate()[0].decode("utf-8").strip()
        # See if chrom contains chr
        chrstart = vcf_chrom.startswith('chr')
        bed_chrom = str(bed_df.iloc[0,0])
        bed_note = bed_chrom.startswith('chr')
        
        if bed_note != chrstart:
            logging.error(f'Warning: Chromosome notations differ between BED file ({bed_chrom}) and VCF/BCF ({vcf_chrom}).')
        # removes or adds "chr" based on analyzed VCF
        bed_df['chr'] = [ norm_chr(chrom, chrstart) for chrom in bed_df['chr'].tolist() ]

        # write to temp file
        bed_df.to_csv(f'{out}_temp.bed', index=False, sep='\t', header=False)

        # gets genotypes at locus of interest, excluding those where 1+ samples missing a genotype call
        bcl_v=f"bcftools view -g ^miss -R {out}_temp.bed {args['<vcf_file>']}"
        
        # Pipe for bcftools
        bcl_view = subprocess.Popen(bcl_v,shell=True, stdout=subprocess.PIPE)
        # splits multiallelic sites into multiple lines
        bcl_norm = subprocess.Popen("bcftools norm -m -",shell=True, stdin=bcl_view.stdout, stdout=subprocess.PIPE)
        bcl_query = subprocess.Popen("bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\n'",shell=True,
         stdin=bcl_norm.stdout, stdout=subprocess.PIPE)
        bcl_query.wait() # Don't do anything else untill bcl_query is done running.

        # output  
        raw_dat = pd.read_csv(StringIO(bcl_query.communicate()[0].decode("utf-8")), sep='\t', 
            header=None, names=['chrom','pos','ref','alt'])
        raw_dat.columns = ['chrom','pos','ref','alt']

        os.remove(f'{out}_temp.bed')
    elif args['--chrom']:
        chrom = args['<locus>']
        print(f'Running on chromosome {chrom}')
        logging.info('Running on entire chromosome. This might take awhile.')
        # get locus info
        # check whether chromosome in VCF file includes "chr" in chromosome
        vcf_chrom = subprocess.Popen(f'bcftools view -H {vcf_in} | cut -f1 | head -1', shell=True, 
            stdout=subprocess.PIPE).communicate()[0].decode("utf-8").strip()

        # See if chrom contains chr
        chrstart = vcf_chrom.startswith('chr')

        chrom = norm_chr(chrom, chrstart)

        bcl_v=f"bcftools view -g ^miss -r {chrom} {args['<vcf_file>']}"
        
        # Pipe for bcftools
        bcl_view = subprocess.Popen(bcl_v,shell=True, stdout=subprocess.PIPE)
        # splits multiallelic sites into multiple lines
        bcl_norm = subprocess.Popen("bcftools norm -m -",shell=True, stdin=bcl_view.stdout, stdout=subprocess.PIPE)
        bcl_query = subprocess.Popen("bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\n'",shell=True,
         stdin=bcl_norm.stdout, stdout=subprocess.PIPE)
        bcl_query.wait() # Don't do anything else untill bcl_query is done running.

        # output  
        raw_dat = pd.read_csv(StringIO(bcl_query.communicate()[0].decode("utf-8")), sep='\t')
        raw_dat.columns = ['chrom','pos','ref','alt']
    else:
        logging.info('Running single locus')

        # get locus info
        locus = args['<locus>']
        chrom = locus.split(':')[0]
        # get locus info
        # check whether chromosome in VCF file includes "chr" in chromosome
        vcf_chrom = subprocess.Popen(f'bcftools view -H {vcf_in} | cut -f1 | head -1', shell=True, 
            stdout=subprocess.PIPE).communicate()[0].decode("utf-8").strip()

        # See if chrom contains chr
        chrstart = vcf_chrom.startswith('chr')

        chrom = norm_chr(chrom, chrstart)

        # properly formatted locus string
        locus=f'{chrom}:'+locus.split(':')[1]

        bcl_v=f"bcftools view -r {locus} {args['<vcf_file>']}"
        
        # Pipe for bcftools
        bcl_view = subprocess.Popen(bcl_v,shell=True, stdout=subprocess.PIPE)
        # splits multiallelic sites into multiple lines
        bcl_norm = subprocess.Popen("bcftools norm -m -",shell=True, stdin=bcl_view.stdout, stdout=subprocess.PIPE)
        bcl_query = subprocess.Popen("bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\n'",shell=True,
         stdin=bcl_norm.stdout, stdout=subprocess.PIPE)
        bcl_query.wait() # Don't do anything else until bcl_query is done running.

        # output  
        raw_dat = pd.read_csv(StringIO(bcl_query.communicate()[0].decode("utf-8")), sep='\t', 
            header=None, names=['chrom','pos','ref','alt'])

    # save to HDF
    out_fname=f'{out}.h5'

    # raw_dat.to_csv(f'{out}.csv')
    raw_dat.to_hdf(out_fname,'all', data_columns=True)

    add_metadata(out_fname, args, os.path.basename(__file__), __version__, 'Gens')

    logging.info('Finished.')


if __name__ == '__main__':
    arguments = docopt(__doc__, version=__version__)
    if arguments['-v']:
        logging.basicConfig(level=logging.INFO, format='[%(asctime)s %(name)s:%(levelname)s ]%(message)s')
    else:
        logging.basicConfig(level=logging.ERROR, format='[%(asctime)s %(name)s:%(levelname)s ]%(message)s')
    print('starting')
    main(arguments)
