#!/bin/bash
#$ -S /bin/bash
#$ -o /netapp/home/kathleen/1kgp_reformat_merge/logs/
#$ -e /netapp/home/kathleen/1kgp_reformat_merge/logs/
#$ -j y
#$ -cwd
#$ -l mem_free=4G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=336:00:00
#$ -t 1-24
#$ -tc 500

export TMPDIR=/scratch
chroms=( null {1..22} X Y )
export chrom=${chroms[$SGE_TASK_ID]}

echo "Running on chromosome ${chrom}."
python /pollard/home/kathleen/projects/AlleleAnalyzer/manuscript_analyses/1000genomes_analysis/src/hg19_analysis/merge_reformatted_variants_1kgp_hg19.py ${chrom}

