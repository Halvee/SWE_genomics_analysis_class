#!/bin/bash

# go to source directory and clone 1kg data to it
mkdir -p src/
cd src/
git clone https://github.com/Halvee/1kg_genotypes_lightweight.git 
cd ../

# make output dir 
mkdir -p results/ancestry_prediction/

# align 1000 genomes phase 3 data with test genotype data
src/1kg_genotypes_lightweight/merge_1000_genomes_with_test_dataset.sh \
data/hapgen_sample3a \
src/1kg_genotypes_lightweight/plink_bedbimfam/g1k_phase3.auto_MAFgt05_genotype_arrays \
bin/plink \
results/ancestry_prediction/test_ref_gts

# get LD-pruned set of SNPs
bin/plink \
--bfile results/ancestry_prediction/test_ref_gts \
--indep 50 5 2 \
--maf 0.05 --mind 0.01 --geno 0.01 --hwe 0.01 \
--out results/ancestry_prediction/test_ref_gts.snps

# do PCA
bin/plink \
--bfile results/ancestry_prediction/test_ref_gts \
--extract results/ancestry_prediction/test_ref_gts.snps.prune.in \
--pca \
--out results/ancestry_prediction/test_ref_gts.pca

# plot PCs
Rscript scripts/plot_plink_eigenvec.R \
results/ancestry_prediction/test_ref_gts.fam \
results/ancestry_prediction/test_ref_gts.pca.eigenvec \
20 \
results/ancestry_prediction/test_ref_gts.pca.pdf 

# run ancestry classification Rscript
Rscript 03_ancestry_prediction.R

exit
