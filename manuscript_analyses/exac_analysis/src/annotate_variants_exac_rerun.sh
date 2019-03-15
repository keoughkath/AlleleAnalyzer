#!/bin/bash
#$ -S /bin/bash
#$ -o /netapp/home/kathleen/exac_annots/logs_rerun/
#$ -e /netapp/home/kathleen/exac_annots/logs_rerun/
#$ -j y
#$ -cwd
#$ -l mem_free=4G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=336:00:00
#$ -t 1-70
#$ -tc 500

export TMPDIR=/scratch
readarray files < /pollard/home/kathleen/projects/AlleleAnalyzer/manuscript_analyses/exac_analysis/dat/non_completed_fnames.txt
readarray regions < /pollard/home/kathleen/projects/AlleleAnalyzer/manuscript_analyses/exac_analysis/dat/regions_nc.txt
regions_to_eval=(null ${regions[@]})
files_to_eval=(null ${files[@]})
echo $SGE_TASK_ID
export file=${files_to_eval[$SGE_TASK_ID]}
echo $file
export region=${regions_to_eval[$SGE_TASK_ID]}
python /pollard/home/kathleen/projects/AlleleAnalyzer/preprocessing/annotate_variants/annot_variants.py -v $file cpf1,SpCas9,SpCas9_VRER,SpCas9_EQR,SpCas9_VQR_1,SpCas9_VQR_2,StCas9,StCas9_2,SaCas9,SaCas9_KKH,nmCas9,cjCas9 /pollard/data/projects/AlleleAnalyzer_data/pam_sites_hg19/ /pollard/data/vertebrate_genomes/human/hg19/hg19/hg19.fa /pollard/data/projects/AlleleAnalyzer_data/exac_data/exac_annotated_variants/$region