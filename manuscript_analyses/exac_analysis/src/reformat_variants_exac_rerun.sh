#!/bin/bash
#$ -S /bin/bash
#$ -o /netapp/home/kathleen/exac_gens/logs3/
#$ -e /netapp/home/kathleen/exac_gens/logs3d/
#$ -j y
#$ -cwd
#$ -l mem_free=4G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=336:00:00
#$ -t 1-69
#$ -tc 500

export TMPDIR=/scratch
echo SGE_TASK_ID
export locus=`sed "${SGE_TASK_ID}q;d" /pollard/home/kathleen/projects/AlleleAnalyzer/manuscript_analyses/exac_analysis/dat/exac_non_comp_regions_for_gens.bed`
echo $locus
chrom=`echo $locus | cut -f1 -d" "`
start=`echo $locus | cut -f2 -d" "`
stop=`echo $locus | cut -f3 -d" "`
name=`echo $locus | cut -f4 -d" "`
export locus_formatted=${chrom}:${start}-${stop}
echo ${locus_formatted}
python /pollard/home/kathleen/projects/AlleleAnalyzer/preprocessing/generate_gens_dfs/get_gens_df.py -v /pollard/data/genetics/ExAc/ExAC.r1.sites.vep.vcf.gz ${locus_formatted} /pollard/data/projects/AlleleAnalyzer_data/exac_data/exac_formatted_variants/$name