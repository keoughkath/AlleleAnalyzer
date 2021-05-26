# AlleleAnalyzer

The CRISPR/Cas system is a highly specific genome editing tool capable of distinguishing alleles differing by even a single base pair. Target sites might carry genetic variations that are not distinguishable by sgRNA designing tools based on one reference genome. AlleleAnalyzer is open-source software that incorporates single nucleotide variants and short insertions and deletions to design sgRNAs for precisely editing one or multiple haplotypes of a sequenced genome *from any species*, currently supporting eleven Cas proteins.  It also leverages patterns of shared genetic variation to optimize sgRNA design for different populations. AlleleAnalyzer is available at https://github.com/keoughkath/AlleleAnalyzer.

**Note**:  If you use this tool and find any issues or anything that is difficult to interpret, please send a message to this Github or to kathleen.keough at gladstone.ucsf.edu

If you use this tool, please cite our [paper in Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1783-3).

## Table of Contents

1. Getting started
* [How it works](https://github.com/keoughkath/ExcisionFinder/wiki/Overview)
* [Install AlleleAnalyzer and its requirements](https://github.com/keoughkath/AlleleAnalyzer/wiki/Install-AlleleAnalyzer-and-its-requirements)
2. Tool descriptions
* [Design allele-specific sgRNAs](https://github.com/keoughkath/AlleleAnalyzer/wiki/Usage:-gen_sgRNAs.py)
* [Using CRISPOR with AlleleAnalyzer](https://github.com/keoughkath/AlleleAnalyzer/wiki/Using-CRISPOR-with-gen_sgRNAs.py)
3. Tutorials
* [Tutorial: Design sgRNAs for allele specific excision of the gene MFN2 in the WTC genome](https://github.com/keoughkath/AlleleAnalyzer/wiki/Tutorial:-Design-sgRNAs-for-allele-specific-excision-of-the-gene-MFN2-in-the-WTC-genome)

## Useful links

[PAM sites in UCSC Genome Browser, hg19](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A61717368-61717468&hgsid=710107349_giMIiVkYz3tMheeUnqbOtAFTIgOo)

[PAM sites in UCSC Genome Browser, hg38](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr11%3A61957117-61957165&hgsid=710108079_SecTcyDrgBPU4AocIPTRF2Uq4Omd)

[UCSC chromosome fastas, hg19.](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/)

[UCSC chromosome fastas, hg38.](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/)


