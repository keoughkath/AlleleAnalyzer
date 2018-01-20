As part of the ExcisionFinder project, we annotated all variants from the 1000 Genomes Project with whether they make, break or are near a PAM site for multiple varieties of Cas, their allele frequency and heterozygosity frequency. If you haven't already, make sure you [download the database](https://github.com/keoughkath/ExcisionFinder/wiki/Download-the-ExcisionFinder-database). 

There are a variety of Python scripts designed to make querying the ExcisionFinder database simple. They are located in the scripts directory of the ExcisionFinder repository that you should have [cloned from Github](https://github.com/keoughkath/ExcisionFinder/wiki/Install-ExcisionFinder-and-its-requirements). Each of these displays use instructions if you use `python script_name.py -h`.

### Query annotated variants for a gene or chromosomal range

getAnnotVariants.py

### Determine targetability for a particular gene based on 1000 Genomes data

getGeneTarg.py

