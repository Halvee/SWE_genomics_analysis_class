

Main <- function(){

  # define input and output files
  in.plink.eigenval <- "results/ancestry_prediction/test_ref_gts.pca.eigenval"
  in.plink.eigenvec <- "results/ancestry_prediction/test_ref_gts.pca.eigenvec"
  refpanel.file <- "src/1kg_genotypes_lightweight/plink_bedbimfam/integrated_call_samples_v3.20130502.ALL.panel"
  outroot <- "results/ancestry_prediction/test_ref_gts.pca.eigenvec.res"
  
  # read 1000 genomes phase 3 reference panel data
  refpanel <- read.table(refpanel.file,
			 header=TRUE, stringsAsFactors=FALSE)
  
  # get european-ancestry reference samples
  eur_refpanel_iids <- subset(refpanel, super_pop == "EUR")$sample
  
  # get eigenvalues from pca and plot
  eval <- scan(in.plink.eigenval, what=numeric())
  pdf(paste0(outroot, ".eigenval.pdf"))
  plot(eval)
  dev.off()
    
  # read eigenvectors from pca and merge with ref panel info
  evec <- read.table(in.plink.eigenvec, header=FALSE, stringsAsFactors=FALSE)
  colnames(evec) <- c("FID","IID", paste0("PC", 1:20))
  evec <- merge(evec, refpanel, by.x="IID", by.y="sample", all.x=TRUE)
  
  # get case/control eigenvectors
  evec.caco <- subset(evec, (IID %in% refpanel$sample) == FALSE)
  evec.caco.pre_pruning <- evec.caco
  evec.ref <- subset(evec, (IID %in% refpanel$sample) == TRUE)
  evec.eur_ref <- subset(evec, IID %in% eur_refpanel_iids)

  # approach 1 : subset on samples that are likely of a single ancestry
  # by defining centroid of known-ancestry samples in reference and removing
  # outliers
  pdf(paste0(outroot, ".pcs_dist.pdf"))
  for (i in 1:4) {
    pc_i <- paste0("PC",i)
    pc_i_eur_ref_mean <- mean(evec.eur_ref[[pc_i]])
    pc_i_eur_ref_sd <- sd(evec.eur_ref[[pc_i]])
    pc_i_range <- c(pc_i_eur_ref_mean - 3*pc_i_eur_ref_sd,
                    pc_i_eur_ref_mean + 3*pc_i_eur_ref_sd)
    evec.caco <- subset(evec.caco,
                        (evec.caco[[pc_i]] >= pc_i_range[1]) &
                        (evec.caco[[pc_i]] <= pc_i_range[2])
                       )
    plot(density(evec.caco.pre_pruning[[pc_i]]),
         col='red',
         main=pc_i,
         xlim=c(min(evec[[pc_i]]), max(evec[[pc_i]])))
    lines(density(evec.eur_ref[[pc_i]]), col='blue')
    abline(v=pc_i_range[1])
    abline(v=pc_i_range[2])
  }
  dev.off()
  cat("N samples (post-EUR subsetting):",nrow(evec.caco), "\n")
  
  # approach 2 : multinomial regression
  # for each sample, return prob of each one having each ancestry label
  # https://bookdown.org/chua/ber642_advanced_regression/multinomial-logistic-regression.html
  # 11.7.3
  evec.ref$super_pop <- factor(evec.ref$super_pop, levels=unique(evec.ref$super_pop))
  library(nnet)
  multi_mo <- multinom(super_pop ~ PC1 + PC2 + PC3 + PC4,
                       data = evec.ref, model = TRUE)
  res <- predict(multi_mo, type='probs',
                 newdata=evec.caco.pre_pruning[,c("IID",paste0("PC",1:4))])
  res <- cbind(data.frame(IID=evec.caco.pre_pruning$IID), res)
  
  # write probabilities to file
  write.table(res,
              file=paste0(outroot, ".ancestry_probs.tsv"),
              row.names=F, col.names=T, sep="\t", quote=F)
  
  # make pdf of pEUR vs pEAS
  pdf(paste0(outroot, ".scatter.pEUR_pEAS.pdf")) 
  plot(jitter(res$EUR, factor=100), jitter(res$EAS, factor=100),
       pch=3, xlab='pEUR', ylab='pEAS') 
  dev.off()
  
  # other approaches possible : 
  # support vector machine (see software 'peddy')
  # random forest model
  
}

if (interactive() == F ){
  Main()
}
