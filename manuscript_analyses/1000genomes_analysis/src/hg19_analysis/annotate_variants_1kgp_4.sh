#!/bin/bash
#$ -S /bin/bash
#$ -o /netapp/home/kathleen/1kgp_hg19_annots/logs/
#$ -e /netapp/home/kathleen/1kgp_hg19_annots/logs/
#$ -j y
#$ -cwd
#$ -l mem_free=4G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=336:00:00
#$ -t 300000-313762
#$ -tc 500

export TMPDIR=/scratch
echo $SGE_TASK_ID
python /pollard/home/kathleen/projects/AlleleAnalyzer/preprocessing/annotate_variants/annot_variants.py -v /pollard/data/projects/AlleleAnalyzer_data/1kgp_data/1kgp_formatted_variants/region_${SGE_TASK_ID}.h5 cpf1,SpCas9,SpCas9_VRER,SpCas9_EQR,SpCas9_VQR_1,SpCas9_VQR_2,StCas9,StCas9_2,SaCas9,SaCas9_KKH,nmCas9,cjCas9 /pollard/data/projects/AlleleAnalyzer_data/pam_sites_hg19/ /pollard/data/vertebrate_genomes/human/hg19/hg19/hg19.fa /pollard/data/projects/AlleleAnalyzer_data/1kgp_data/1kgp_annotated_variants/region${SGE_TASK_ID}