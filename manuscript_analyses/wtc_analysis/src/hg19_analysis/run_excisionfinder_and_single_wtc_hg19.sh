#!/bin/bash
#$ -S /bin/bash
#$ -o /netapp/home/kathleen/excisionfinder_wtc_hg19_s/
#$ -e /netapp/home/kathleen/excisionfinder_wtc_hg19_s/
#$ -j y
#$ -cwd
#$ -l mem_free=6G
#$ -l scratch=1G
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

python /pollard/home/kathleen/projects/AlleleAnalyzer/scripts/ExcisionFinder.py -v /pollard/home/kathleen/projects/AlleleAnalyzer/manuscript_analyses/1000genomes_analysis/get_gene_list/gene_list_hg19.tsv ${gene} /pollard/data/projects/AlleleAnalyzer_data/wtc_data/hg19/wtc_annotated_variants_by_chrom2/${chrom}.h5 10000 cpf1,SpCas9,SpCas9_VRER,SpCas9_EQR,SpCas9_VQR_1,SpCas9_VQR_2,StCas9,StCas9_2,SaCas9,SaCas9_KKH,nmCas9,cjCas9 /pollard/home/kathleen/conklin_wt_seq_data/wtc_wgs_data/phased_yin/wtc_PASS_hg19.phased.vcf.gz /pollard/data/projects/AlleleAnalyzer_data/wtc_data/hg19/wtc_excisionfinder_results2/${chrom}_ef_results/${gene}

python /pollard/home/kathleen/projects/AlleleAnalyzer/manuscript_analyses/wtc_analysis/src/hg19_analysis/single_cut_targ_wtc_hg19.py -v /pollard/home/kathleen/projects/AlleleAnalyzer/manuscript_analyses/1000genomes_analysis/get_gene_list/gene_list_hg19.tsv ${gene} /pollard/data/projects/AlleleAnalyzer_data/wtc_data/hg19/wtc_annotated_variants_by_chrom2/${chrom}.h5 cpf1,SpCas9,SpCas9_VRER,SpCas9_EQR,SpCas9_VQR_1,SpCas9_VQR_2,StCas9,StCas9_2,SaCas9,SaCas9_KKH,nmCas9,cjCas9 /pollard/home/kathleen/conklin_wt_seq_data/wtc_wgs_data/phased_yin/wtc_PASS_hg19.phased.vcf.gz /pollard/data/projects/AlleleAnalyzer_data/wtc_data/hg19/wtc_single_targ2/

qstat -j $JOB_ID
