With ExcisionFinder, we've provided the processed analyses of the 1000 Genomes cohort. Here are the instructions for how to apply ExcisionFinder to your own genome(s) of interest to identify allele-specific CRISPR sites. 

Before you run ExcisionFinder, you need to generate some input files, (1) variant targetability and (2) genotypes dataframes. For this, you need a phased VCF or BCF file. You can find more information on phasing [here](https://faculty.washington.edu/browning/beagle/beagle.html). 

### Generate explicit genotype files

Just plug your phased VCF/BCF file into get_chr_tables.sh, doing each chromosome separately, even if your VCF/BCF has all the chromosomes together. This script will figure it out. Note that this bash script depends on its companion python script, fix_chr_tables.py, which is in the same directory. The Python script is essential to properly splitting the multiallelic sites. 

*Sample Usage*:

bash get_chr_tables.sh wtc_chr1_phased.vcf.gz WTC 1 wtc_chr1_test/

In the above sample run command, wtc_chr1_phased.vcf.gz is the phased VCF file for our sample genome of interest, WTC. Please refer to the ExcisionFinder manuscript for more details on WTC. WTC is the sample name. If you have a samples file, you will have to refer to the [bcftools documentation](http://www.htslib.org/doc/bcftools.html) to make a small change to the bcftools command used in the bash script. Finally, wtc_chr1_test/ here is the output directory. 

### Generate variant targetability files

This part annotated variants that make, break or are near PAM sites (variant targetability). This part requires that you have downloaded the pre-computed locations of PAM sites for all varieties of Cas analyzed by ExcisionFinder (available on the [figshare page](https://figshare.com/account/home#/projects/26080) for this project). Additionally, you will need the fasta file for GrCh37 (hg19), also available on the [figshare page](https://figshare.com/account/home#/projects/26080). 

*Sample Usage*:

python gen_targ_dfs.py chr1_gens.hdf5 1 hg19_pams/ grch37.fa wtc_chr1_test/

In the above sample run command, chr1_gens.hdf5 is the chr1 gens file for WTC that was generated using the above command with get_chr_tables.sh. 1 is the chromosome, hg19_pams/ is the directory containing PAM site locations in hg19, grch37.fa is the GrCh37 (hg19) fasta file, and wtc_chr1_test/ is the output file where the targ_dfs will go.


### Analyzing gene targetability with ExcisionFinder

This is the part that annotates targetability for each gene in each individual in your cohort, and additionally puts out dataframes annotating targetability of each haplotype for each individual at each variant position. Note: it requires that the custom Python module crisprtools (available on this github) is in the current working directory.

*Sample Usage*:

ExcisionFinder.py -v gene_annots_wsize DEFB128 20 targ_dfs/ gen_dfs/ outdir/

In the above sample, gene_annots_wsize is the gene annotations file used for the EF manuscript, available on [figshare](https://figshare.com/account/home#/projects/26080). You can also make your own if you're interested in a particular set of gene, just ensure it follows the same format. DEFB128 is my favorite test gene, as it is small and runs relatively quickly, but this is where whichever gene being analyzed is put. A small bash script can parallelize this process and analyze multiple genes concurrently. 20 indicates the chromosome on which the gene resides. targ_dfs/ is the location of the targetability dfs from gen_targ_dfs.py and outdir/ is where the output file, i.e. in this case chr20_out.hdf5, will be saved. This particular command will output a file chr20_out.hdf5 in the output directory as ExcisionFinder is designed to put all genes from a particular chromosome together. 