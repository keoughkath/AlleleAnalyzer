#!/bin/bash
#$ -S /bin/bash
#$ -o /netapp/home/kathleen/wtc_hg19_annots_qb3/logs2/
#$ -e /netapp/home/kathleen/wtc_hg19_annots_qb3/logs2/
#$ -j y
#$ -cwd
#$ -l mem_free=4G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=336:00:00
#$ -t 60000-100000
#$ -tc 500

export TMPDIR=/scratch
echo $SGE_TASK_ID
python /pollard/home/kathleen/projects/AlleleAnalyzer/preprocessing/annotate_variants/annot_variants.py -v /pollard/data/projects/AlleleAnalyzer_data/wtc_data/hg19/wtc_formatted_variants/region_${SGE_TASK_ID}.h5 cpf1,SpCas9,SpCas9_VRER,SpCas9_EQR,SpCas9_VQR_1,SpCas9_VQR_2,StCas9,StCas9_2,SaCas9,SaCas9_KKH,nmCas9,cjCas9 /pollard/data/projects/AlleleAnalyzer_data/pam_sites_hg19/ /pollard/data/vertebrate_genomes/human/hg19/hg19/hg19.fa /netapp/home/kathleen/wtc_hg19/wtc_annotated_variants/region_${SGE_TASK_ID}