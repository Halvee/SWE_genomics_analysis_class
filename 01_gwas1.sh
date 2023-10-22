#!/bin/bash

# data paths
BEDBIMFAM=data/sim_sim1a_eur_sa_merge.miss
OUTDIR=results/gwas1
OUTROOT=$OUTDIR/sim_sim1a_eur_sa_merge

# we will follow the steps that are used by default by Ricopili pipeline        
# (https://sites.google.com/a/broadinstitute.org/ricopili/preimputation-qc),    
# the analysis pipeline that is used by the Psychiatric Genomics Consortium 

# first, make your output directory
mkdir -p $OUTDIR

# what happens if we try to run a case/control assoc test right now without
# any quality control?
bin/plink \
--bfile $BEDBIMFAM \
--assoc \
--out $OUTROOT.caco
# and then make a QQ plot
Rscript scripts/qq.plink.R $OUTROOT.caco.assoc $OUTROOT.caco.assoc.png

# first let's perform missingness-based QC. 
# ID variants with high missingness (>0.05), exclude, and then find 
# persons where even excluding these variants missingness is high (>0.02)
bin/plink \
--bfile $BEDBIMFAM \
--missing \
--out $OUTROOT

# get snp ids with missingness > 0.05
awk '{if ($5 > 0.05) {print $2}}' \
$OUTROOT.lmiss \
> $OUTROOT.lmiss_gt_05.snpid.list

# exclude samples with high rates of missingness (>0.02) along the snps that 
# don't have high missingness
bin/plink \
--bfile $BEDBIMFAM \
--exclude $OUTROOT.lmiss_gt_05.snpid.list \
--mind 0.02 \
--make-bed \
--out $OUTROOT.miss1_pruned

# compare missingness rates in cases versus controls
bin/plink \
--bfile $OUTROOT.miss1_pruned \
--test-missing \
--out $OUTROOT.miss1_pruned.caco_missing

# get SNPs where case/control missingness rate doesn't differ by more than 
# 0.02 or greater
awk '{diff=$4-$3; if ((-0.02 > diff) || (diff > 0.02)) {print $2}}' \
$OUTROOT.miss1_pruned.caco_missing.missing \
> $OUTROOT.miss1_pruned.caco_missing_gt_02.snpid.list

# subset  data on SNPs that meet the following criteria:
# 1. missingness across full cohort <= 0.02
# 2. case vs control difference in missingness <= 0.02
bin/plink \
--bfile $OUTROOT.miss1_pruned \
--exclude $OUTROOT.miss1_pruned.caco_missing_gt_02.snpid.list \
--mind 0.02 \
--make-bed \
--out $OUTROOT.miss2_pruned

# run plink function for hardy-weinberg equilibrium
bin/plink \
--bfile $OUTROOT.miss2_pruned \
--hardy \
--out $OUTROOT.miss2_pruned.hardy

# get SNPs that fail HWE in cases (10^-10)
python scripts/pval_to_log10p.py \
  $OUTROOT.miss2_pruned.hardy.hwe \
  9 \
| awk '{if (($3=="AFF") && ($9 >10)) {print $2}}' \
> $OUTROOT.miss2_pruned.hardy.ca_fail.snpid.list

# get SNPs that fail HWE in controls (10^-6)
python scripts/pval_to_log10p.py \
  $OUTROOT.miss2_pruned.hardy.hwe \
  9 \
| awk '{if (($3=="UNAFF") && ($9 >6)) {print $2}}' \
> $OUTROOT.miss2_pruned.hardy.co_fail.snpid.list

# combine case and control HWE fail snp sets
cat $OUTROOT.miss2_pruned.hardy.ca_fail.snpid.list \
    $OUTROOT.miss2_pruned.hardy.co_fail.snpid.list \
| sort \
| uniq \
> $OUTROOT.miss2_pruned.hardy.caco_fail.snpid.list

# remove snps that fail HWE
bin/plink \
--bfile $OUTROOT.miss2_pruned \
--exclude $OUTROOT.miss2_pruned.hardy.caco_fail.snpid.list \
--make-bed \
--out $OUTROOT.miss2_hwe_pruned

# what happens if we try to run a case/control assoc test right now without
# any quality control?
bin/plink \
--bfile $OUTROOT.miss2_hwe_pruned \
--assoc \
--out $OUTROOT.miss2_hwe_pruned.caco
# and then make a QQ plot
Rscript scripts/qq.plink.R $OUTROOT.miss2_hwe_pruned.caco.assoc \
$OUTROOT.miss2_hwe_pruned.caco.assoc.png
