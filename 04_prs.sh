#!/bin/bash

# make output dir
mkdir -p results/prs/Height_PRS/

DL=1
if [[ $DL == 1 ]]
then

# go to output dir
cd results/prs/

# get GWAS summary statistics for height (Yengo 2022)
wget -N \
"https://504394d8-624a-4827-9f25-95a83cd9675a.filesusr.com/archives/1d0101_394c2d3120ba4d0a9f6326ed56ff8854.gz?dn=GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz"
cp \
1d0101_394c2d3120ba4d0a9f6326ed56ff8854.gz\?dn\=GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz \
GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz

# get pre-computed score weights for height (Yengo 2022)
wget -N \
"https://504394d8-624a-4827-9f25-95a83cd9675a.filesusr.com/archives/1d0101_ea67ab206da4462f9cbca90a6b4c011f.gz?dn=GIANT_HEIGHT_YENGO_2022_PGS_WEIGHTS_EUR.gz"
cp \
1d0101_ea67ab206da4462f9cbca90a6b4c011f.gz\?dn\=GIANT_HEIGHT_YENGO_2022_PGS_WEIGHTS_EUR.gz \
GIANT_HEIGHT_YENGO_2022_PGS_WEIGHTS_EUR.gz

# get ld information for EUR, preformatted for PRScs
# (we use 1000 genomes here, but UK biobank LD is also available, will be more accurate)
wget -N \
"https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz"
tar -xvzf ldblk_1kg_eur.tar.gz ldblk_1kg_eur/ldblk_1kg_chr22.hdf5
tar -xvzf ldblk_1kg_eur.tar.gz ldblk_1kg_eur/snpinfo_1kg_hm3

# back to parent directory
cd ../../

fi

# format sumstats for input into PRS calculation
# plink_CT : need SNP, P-value columns
# PRS-CS : (see https://github.com/getian107/PRScs#using-prs-cs)
echo -e "SNP\tA1\tA2\tBETA\tSE\tP" \
> results/prs/Height.Yengo_2022.sum_stats_file.tsv
# input table columns :
# SNPID	RSID	CHR	POS	EFFECT_ALLELE	OTHER_ALLELE	EFFECT_ALLELE_FREQ	BETA	SE	P	N
gunzip -c results/prs/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz \
| awk '{OFS="\t"; if (NR!=1) {print $2,$5,$6,$8,$9,$10}}' \
>>results/prs/Height.Yengo_2022.sum_stats_file.tsv

# derive subset of SNPs that are nonambiguous (no A/T, G/C)
# Reason : with A/T and G/C SNPs you could difference in what strand a call is relative to
awk '{A1A2["AG"]=1;A1A2["AC"]=1;A1A2["CA"]=1;A1A2["CT"]=1;A1A2["GA"]=1;A1A2["GT"]=1;A1A2["TC"]=1;A1A2["TG"]=1;
      ALLELES=$2$3;
     if ((NR==1) || (ALLELES in A1A2)) {print $0}}' \
results/prs/Height.Yengo_2022.sum_stats_file.tsv \
> results/prs/Height.Yengo_2022.sum_stats_file.nonambig.tsv 

# lastly, let's make a file specific for PRS-CS, which only wants the first 5
# columns in this file
cut -f 1-5 \
results/prs/Height.Yengo_2022.sum_stats_file.nonambig.tsv \
> results/prs/Height.Yengo_2022.sum_stats_file.nonambig.PRS-CS.tsv

RUN_PLINK_CT=1
if [[ $RUN_PLINK_CT == 1 ]]
then

# see https://choishingwan.github.io/PRS-Tutorial/plink/ for more details
# as well as https://zzz.bwh.harvard.edu/plink/clump.shtml

# step 1 : LD clumping : reduce sumstats to independant loci represented by most
# significant SNPs in each LD block
# note 1 : cohort should be post-qc, single-ancestry if using for clumping
# note 2 : if cohort <500, then use ref panel genotypes for clumping instead
bin/plink \
--bfile data/sim_sim2a_eur_sa_merge.miss \
--chr 22 \
--clump-p1 0.05 \
--clump-r2 0.1 \
--clump-kb 250 \
--clump results/prs/Height.Yengo_2022.sum_stats_file.nonambig.tsv  \
--clump-snp-field SNP \
--clump-field P \
--out results/prs/Height_PRS/Height_PRS.plink_CT

# get index SNPs from LD clumps 
awk '{print $3}' \
results/prs/Height_PRS/Height_PRS.plink_CT.clumped \
> results/prs/Height_PRS/Height_PRS.plink_CT.clumped.snpid.list

fi

RUN_PRSCS=1
if [[ $RUN_PRSCS == 1 ]]
then

# run PRScs command to get fitted weights (chr22 only for sake of time)
# n_gwas = 4,080,687 (EUR subset)
python src/PRScs/PRScs.py \
--sst_file results/prs/Height.Yengo_2022.sum_stats_file.nonambig.PRS-CS.tsv \
--n_gwas 4080687 \
--ref_dir results/prs/ldblk_1kg_eur \
--bim_prefix data/sim_sim2a_eur_sa_merge.miss \
--chrom 22 \
--out_dir results/prs/Height_PRS/Height_PRS

fi

# copy results to new file
cat \
results/prs/Height_PRS/Height_PRS_pst_eff_a1_b0.5_phiauto_chr*.txt \
> results/prs/Height_PRS/Height_PRS_pst_eff_a1_b0.5_phiauto.autosome.txt

# use plink to calculate PRS per sample using computed SNP weights (plink_CT)
# --score <filename> [variant ID col.] [allele col.] [score col.] 
bin/plink \
--bfile data/sim_sim2a_eur_sa_merge.miss \
--extract results/prs/Height_PRS/Height_PRS.plink_CT.clumped.snpid.list \
--score results/prs/Height.Yengo_2022.sum_stats_file.nonambig.PRS-CS.tsv \
        1 2 4 \
--out results/prs/Height_PRS/sim_sim2a_eur_sa_merge.miss.height_prs.plink_CT

# use plink to calculate PRS per sample using computed SNP weights (PRS-CS)
bin/plink \
--bfile data/sim_sim2a_eur_sa_merge.miss \
--score results/prs/Height_PRS/Height_PRS_pst_eff_a1_b0.5_phiauto.autosome.txt \
        2 4 6 \
--out results/prs/Height_PRS/sim_sim2a_eur_sa_merge.miss.height_prs.PRS-CS 

# LDpred2 : another program for calculating PRS without p-value thresholding
# https://privefl.github.io/bigsnpr/articles/LDpred2.html
# Tutorial video : https://www.youtube.com/watch?v=aya8WsNAu6U
# LDpred2-auto

exit
