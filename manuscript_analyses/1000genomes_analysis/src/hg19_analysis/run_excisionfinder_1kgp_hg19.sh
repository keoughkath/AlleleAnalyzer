#!/bin/bash
#$ -S /bin/bash
#$ -o /netapp/home/kathleen/excisionfinder_1kgp_hg19/logs/
#$ -e /netapp/home/kathleen/excisionfinder_1kgp_hg19/logs/
#$ -j y
#$ -cwd
#$ -l mem_free=4G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=336:00:00
#$ -t 1-27091
#$ -tc 500

export TMPDIR=/scratch
readarray all_genes < /pollard/home/kathleen/projects/AlleleAnalyzer/manuscript_analyses/1000genomes_analysis/get_gene_list/genes_hg19.txt
genes=(null ${all_genes[@]})
export gene=${genes[$SGE_TASK_ID]}

readarray all_chroms < /pollard/home/kathleen/projects/AlleleAnalyzer/manuscript_analyses/1000genomes_analysis/get_gene_list/chroms_hg19.txt
chroms=(null ${all_chroms[@]})
export chrom=${chroms[$SGE_TASK_ID]}

echo $SGE_TASK_ID
echo ${gene}
echo ${chrom}

python /pollard/home/kathleen/projects/AlleleAnalyzer/scripts/ExcisionFinder.py -v /pollard/home/kathleen/projects/AlleleAnalyzer/manuscript_analyses/1000genomes_analysis/get_gene_list/gene_list_hg19.tsv ${gene} /pollard/data/projects/AlleleAnalyzer_data/1kgp_data/hg19_analysis/1kgp_annotated_variants_by_chrom/${chrom}_annotated.h5 10000 cpf1,SpCas9,SpCas9_VRER,SpCas9_EQR,SpCas9_VQR_1,SpCas9_VQR_2,StCas9,StCas9_2,SaCas9,SaCas9_KKH,nmCas9,cjCas9 /pollard/data/genetics/1kg/phase3/hg19/ALL.${chrom}.*.bcf /pollard/data/projects/AlleleAnalyzer_data/1kgp_data/hg19_analysis/1kgp_excisionfinder_results/results_by_chrom/${chrom}_ef_results/${gene}

qstat -j $JOB_ID
