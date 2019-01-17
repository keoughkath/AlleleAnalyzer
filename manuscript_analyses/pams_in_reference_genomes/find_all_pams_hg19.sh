#!/bin/bash
#$ -S /bin/bash
#$ -o /netapp/home/kathleen/hg19_pams/logs/
#$ -e /netapp/home/kathleen/hg19_pams/logs/
#$ -j y
#$ -cwd
#$ -l mem_free=4G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=336:00:00
#$ -t 1-24
#$ -tc 500

export TMPDIR=/scratch
export chroms=( null {1..22} Y X )
export chrom=${chroms[$SGE_TASK_ID]}
echo $SGE_TASK_ID
echo ${chrom}

time python /pollard/home/kathleen/projects/AlleleAnalyzer/preprocessing/find_pams_in_reference/pam_pos_genome.py chr${chrom} /pollard/data/vertebrate_genomes/human/hg19/hg19/hg19.fa cpf1,SpCas9,SpCas9_VRER,SpCas9_EQR,SpCas9_VQR_1,SpCas9_VQR_2,StCas9,StCas9_2,SaCas9,SaCas9_KKH,nmCas9,cjCas9 /pollard/data/projects/AlleleAnalyzer_data/pam_sites_hg19/
qstat -j $JOB_ID