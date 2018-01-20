ExcisionFinder is a database and toolset for analysis and design of allele-specific CRISPR editing that leverages >2,500 human genomes to estimate feasibility of allele-specific loci in varied individuals and design allele-specific sgRNAs in individual genomes. 

Given a gene or user-defined region, ExcisionFinder can tell you how many people are allele-specifically targetable at that gene and give you information about which populations are most targetable, identify sgRNA pairs that are shared by large groups of people, and identify combinations of sgRNA pairs that would cover the most people using data from the [1000 Genomes Project](http://www.internationalgenome.org/).

Given a VCF file for a human genome, ExcisionFinder can determine which genes are targetable in that specific individual, and design allele-specific sgRNAs for a gene and/or user-defined region of choice. 

# How it works

![](https://github.com/keoughkath/ExcisionFinder/blob/master/imgs/pipeline.jpg)

**Reference dataset**: ExcisionFinder relies on data from the [1000 Genomes Project](http://www.internationalgenome.org/) to explore targetability of a gene or genomic region in varied people. 

**Analysis of individual genomes**: Given a VCF, ExcisionFinder annotates the variants in the individual (based on human assembly genome hg19) and identifies allele-specific CRISPR sites genome-wide. These can then be used to design allele-specific sgRNAs. 

