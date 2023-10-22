#!/bin/bash

# data paths
BEDBIMFAM=results/gwas1/sim_sim1a_eur_sa_merge.miss2_hwe_pruned
OUTDIR=results/gwas2
OUTROOT=$OUTDIR/sim_sim1a_eur_sa_merge

# make sure output directory exists
mkdir -p $OUTDIR

# check heterozygosity in samples
bin/plink \
--bfile $BEDBIMFAM \
--het \
--out $OUTROOT.miss2_hwe_pruned.heterozygosity

# plot the distribution of heterozygosity along with the proposed threshold
Rscript scripts/values_density_plot.R \
$OUTROOT.miss2_hwe_pruned.heterozygosity.het \
F \
$OUTROOT.miss2_hwe_pruned.heterozygosity.het.density_plot.png \
"neg0.2,0.2"

# flag samples for removal where heterozygosity not between (-0.2, 0.2)
awk '{if (($6 < -0.2) || (0.2 < $6)) {print $1,$2}}' \
$OUTROOT.miss2_hwe_pruned.heterozygosity.het \
> $OUTROOT.miss2_hwe_pruned.heterozygosity.fhet_fail.fidiid.txt

# remove heterozygosity-failing samples from dataset
bin/plink \
--bfile $BEDBIMFAM \
--remove $OUTROOT.miss2_hwe_pruned.heterozygosity.fhet_fail.fidiid.txt \
--make-bed \
--out $OUTROOT.miss2_hwe_het_pruned

# get LD-pruned autosomal set of SNPs
# https://zzz.bwh.harvard.edu/plink/summary.shtml
bin/plink \
--bfile $OUTROOT.miss2_hwe_het_pruned \
--indep 50 5 2 \
--out $OUTROOT.miss2_hwe_het_pruned.snpset

# let's do pairwise relatedness checks on all samples. Any cryptic relatedness?
bin/plink \
--bfile $OUTROOT.miss2_hwe_het_pruned \
--extract $OUTROOT.miss2_hwe_het_pruned.snpset.prune.in \
--genome gz \
--out $OUTROOT.miss2_hwe_het_pruned.relatedness

# any sample pairs with 1st or 2nd degree relatedness evidence? PI_HAT > 0.2
gunzip -c $OUTROOT.miss2_hwe_het_pruned.relatedness.genome.gz \
| awk '{if ((NR!=1) && ($10 > 0.2)) {print $1,$2; print $3,$4}}' \
> $OUTROOT.miss2_hwe_het_pruned.relatedness.genome.pihat_gt_02.fidiid.txt

# no samples with cryptic relatedness based on the output above

# do PCA
bin/plink \
--bfile $OUTROOT.miss2_hwe_het_pruned \
--extract $OUTROOT.miss2_hwe_het_pruned.snpset.prune.in \
--pca \
--out $OUTROOT.miss2_hwe_het_pruned.pca

# plot PCs
Rscript scripts/plot_plink_eigenvec.R \
$OUTROOT.miss2_hwe_het_pruned.fam \
$OUTROOT.miss2_hwe_het_pruned.pca.eigenvec \
4 \
$OUTROOT.miss2_hwe_het_pruned.pca.pdf

# density plot for PC1
Rscript scripts/values_density_plot.R \
$OUTROOT.miss2_hwe_het_pruned.pca.eigenvec \
3 \
$OUTROOT.miss2_hwe_het_pruned.pca.PC1.density_plot.png

# get rid of samples with PC1<0 as these are outliers (likely non-EUR).
awk '{if ($3 < 0) {print $1,$2}}' \
$OUTROOT.miss2_hwe_het_pruned.pca.eigenvec \
> $OUTROOT.miss2_hwe_het_pruned.pca.PC1_lt_0.fidiid.txt

# make bedbimfam file with these samples removed
bin/plink \
--bfile $OUTROOT.miss2_hwe_het_pruned \
--remove $OUTROOT.miss2_hwe_het_pruned.pca.PC1_lt_0.fidiid.txt \
--make-bed \
--out $OUTROOT.miss2_hwe_het_nonEUR_pruned

# redo PCA
bin/plink \
--bfile $OUTROOT.miss2_hwe_het_nonEUR_pruned \
--extract $OUTROOT.miss2_hwe_het_pruned.snpset.prune.in \
--pca \
--out $OUTROOT.miss2_hwe_het_nonEUR_pruned.pca

# plot PCs
Rscript scripts/plot_plink_eigenvec.R \
$OUTROOT.miss2_hwe_het_nonEUR_pruned.fam \
$OUTROOT.miss2_hwe_het_nonEUR_pruned.pca.eigenvec \
4 \
$OUTROOT.miss2_hwe_het_nonEUR_pruned.pca.pdf

# do case/control association with QQ plot
bin/plink \
--bfile $OUTROOT.miss2_hwe_het_nonEUR_pruned  \
--assoc \
--out $OUTROOT.miss2_hwe_het_nonEUR_pruned.caco 
# and then make QQ plot
Rscript scripts/qq.plink.R \
$OUTROOT.miss2_hwe_het_nonEUR_pruned.caco.assoc \
$OUTROOT.miss2_hwe_het_nonEUR_pruned.caco.assoc.png

# get list of genomewide significant snps
awk '{if ((NR!=1) && ($9 < 5e-8)) {print $2}}' \
$OUTROOT.miss2_hwe_het_nonEUR_pruned.caco.assoc \
> $OUTROOT.miss2_hwe_het_nonEUR_pruned.caco.assoc.gws.snpid.list

# get caco missingness and hwe per snp                                          
bin/plink \
--bfile $OUTROOT.miss2_hwe_het_nonEUR_pruned \
--extract $OUTROOT.miss2_hwe_het_nonEUR_pruned.caco.assoc.gws.snpid.list \
--missing \
--test-missing \
--hardy \
--freq \
--assoc \
--out $OUTROOT.miss2_hwe_het_nonEUR_pruned.caco.assoc.gws

# make new dataset with higher missingness stringency
bin/plink \
--bfile $OUTROOT.miss2_hwe_het_nonEUR_pruned  \
--geno 0.01 \
--make-bed \
--out $OUTROOT.miss3_hwe_het_nonEUR_pruned

# do one more association test, this time with higher stringency missingness threshold
bin/plink \
--bfile $OUTROOT.miss3_hwe_het_nonEUR_pruned  \
--assoc \
--out $OUTROOT.miss3_hwe_het_nonEUR_pruned.caco 
# and then make QQ plot
Rscript scripts/qq.plink.R \
$OUTROOT.miss3_hwe_het_nonEUR_pruned.caco.assoc \
$OUTROOT.miss3_hwe_het_nonEUR_pruned.caco.assoc.png
