

Main <- function(){
  ARGS <- commandArgs(trailingOnly=T)
  if (length(ARGS) != 4) {
    cat("plot_plink_eigenvec.R <in.fam> <in.eigenvec> ",
        "<n_pcs_to_plot> <out.pdf>\n")
    q()
  }
  in.fam <- ARGS[1]
  in.eigenvec <- ARGS[2]
  n_pcs_to_plot <- as.numeric(ARGS[3])
  out.pdf <- ARGS[4]

  # read files
  fam <- read.table(in.fam, stringsAsFactors=F)
  colnames(fam) <- c("FID","IID","PAT","MAT","SEX","PHE")
  evec <- read.table(in.eigenvec, stringsAsFactors=F)
  colnames(evec) <- c("FID","IID",paste0("PC",1:20))

  # get case iids
  ca_iids <- subset(fam, PHE==2)$IID

  # open filehandle to output pdf
  pdf(out.pdf)

  # for each set of pcs (i, i+1) ..
  i <- 1
  while (i < n_pcs_to_plot) {

    # derive j as i+1
    j <- i + 1

    # derive pc i and j
    pc_i <- paste0("PC",i)
    pc_j <- paste0("PC",j)

    # plot pc i and j
    plot(evec[[pc_i]], evec[[pc_j]], xlab=pc_i, ylab=pc_j)

    # get case subset of evec and plot
    evec.ca <- subset(evec, IID %in% ca_iids)
    points(evec.ca[[pc_i]], evec.ca[[pc_j]], col='red')

    # increment i
    i <- i + 2
  }

  # close filehandle to output pdf
  dev.off()

}

if (interactive() == F) {
  Main()
}
