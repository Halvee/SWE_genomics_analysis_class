

Main <- function(){
  ARGS <- commandArgs(trailingOnly=T)
  if ((length(ARGS) != 3) && (length(ARGS) != 4)) {
    cat("values_density_plot.R <table.txt> <column_number|column_name> ",
        "<out.png> [value_threshold1,..,value_thresholdN]\n")
    q()
  }
  table.txt<-ARGS[1]
  column<-ARGS[2]
  out.png<-ARGS[3]
  thresholds<-c(NA)
  if (length(ARGS) == 4) {
    thresholds <- ARGS[4]
    thresholds <- gsub("neg","-",thresholds)
    thresholds <- strsplit(thresholds, ",")[[1]]
    thresholds <- as.numeric(thresholds)
  }

  # get values
  df <- read.table(table.txt, header=T, stringsAsFactors=F)
  is_not_numeric <- suppressWarnings(is.na(as.numeric(column)))
  if (is_not_numeric == F) {
    column <- as.numeric(column)
  }
  values <- df[, column]
  
  # make density plot
  png(file=out.png, res=300, height = 2000, width = 2000)
  plot(density(values), main=column)
  if (is.na(thresholds[1]) == F) {
    for (threshold in thresholds) {
      abline(v=threshold, col='red')
    }
  }
  dev.off()
}

if (interactive() == F) {
  Main()
}
