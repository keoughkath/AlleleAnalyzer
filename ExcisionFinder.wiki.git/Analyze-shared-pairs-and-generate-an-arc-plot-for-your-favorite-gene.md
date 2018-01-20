In order to complete this tutorial, you must have cloned the ExcisionFinder repository. 

For this tutorial, we will analyze the gene *BEST1*, a gene that encodes for an important protein in the retina which, when mutated, can cause dominant negative retinal dystrophy.

### 1.) Download the necessary data.

If you're using ExcisionFinder on your laptop, you won't want or need the entire database for all genes in the human genome, so you can just download what you need. *BEST1* is on chromosome 11, so we can go [download the necessary files only for chromosome 11](filepath). For simplicity, create a directory (`mkdir best1_ef_dat`) that you can use to keep it all in one place.

The files you will need to have downloaded are: 

*all_samples.npy
*gene_annots_w_size
*sample_legend.tsv
*hap_targ_ind_11.hdf5 (in gene_targ_outputs/chr11/)
*chr11_targ.hdf5 (in targ_dfs)

### 2.) Use the python script gen_arcplot_df_genewrange.py to generate a .tsv file detailing sharing of allele-specific CRISPR sites.

For this tutorial, we're going to analyze the overall 1000 genomes cohort for all types of cas first, then specify later. 

Navigate to the plotting scripts directory of ExcisionFinder:

`cd ExcisionFinder/paper_figures/plotting_scripts`

Generate the shared pairs file:

`python gen_arcplot_df_genewrange.py BEST1 0 --type=a 10000 