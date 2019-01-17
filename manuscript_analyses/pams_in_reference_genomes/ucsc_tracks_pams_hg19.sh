echo "Making tracks for each cas"

python ucsc_tracks_pams_hg19.py

echo "Sorting the BED file"

parallel -j10 " sort -k1,1 -k2,2n /pollard/data/projects/AlleleAnalyzer_data/pam_beds_hg19/ucsc_browser/{}_ucsc.bed > /pollard/data/projects/AlleleAnalyzer_data/pam_beds_hg19/ucsc_browser/{}_ucsc.sorted.bed " ::: SpCas9 SpCas9_VRER SpCas9_EQR SpCas9_VQR_1 SpCas9_VQR_2 StCas9 StCas9_2 SaCas9 SaCas9_KKH nmCas9 cjCas9 cpf1

echo "Generating bigBed file for each cas"

parallel -j10 " bedToBigBed /pollard/data/projects/AlleleAnalyzer_data/pam_beds_hg19/ucsc_browser/{}_ucsc.sorted.bed /pollard/data/vertebrate_genomes/human/hg19/hg19/hg19.sizes /pollard/data/projects/AlleleAnalyzer_data/pam_beds_hg19/ucsc_browser/{}_ucsc.bb " ::: SpCas9 SpCas9_VRER SpCas9_EQR SpCas9_VQR_1 SpCas9_VQR_2 StCas9 StCas9_2 SaCas9 SaCas9_KKH nmCas9 cjCas9 cpf1