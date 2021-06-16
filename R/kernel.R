#' Kernel
#'
#' Description of function
#' @param data450 An HM450 data set with CpG sites as columns
#' @param data850 An EPIC probe data set with CpG sites as columns
#' @param chrom Indicator of chromosome
#' @param path File path for output file export
#' @param sigma Optionally fix sigma to desired value
#' @export
kernel <- function(data450, data850, chrom, path, sigma = NULL) {
  # special case for sex chromosomes
  chrom_val <- NULL
  if (chrom == 'X') {
    chrom_val <- 23
  }
  if (chrom == 'Y') {
    chrom_val <- 24
  }
  else {
    chrom_val <- chrom
  }

   # impute missing values with the column mean (source data and user data)
  for(i in 1:ncol(data850)) {
    data850[is.na(data850[,i]), i] <- mean(data850[,i], na.rm = TRUE)
  }

  for(i in 1:ncol(data450)) {
    data450[is.na(data450[,i]), i] <- mean(data450[,i], na.rm = TRUE)
  }

  # split data850 into x_train (450k feature set) and y_train (EPIC only feature set)
  indicator <- metImpute::indicatordata
  ind_xtrain <- subset(indicator, Methyl450_Loci == TRUE)
  ind_ytrain <- subset(indicator, Methyl450_Loci == FALSE)
  x_train <- subset(data850, colnames(data850) %in% ind_xtrain$IlmnID)
  y_train <- subset(data850, colnames(data850) %in% ind_ytrain$IlmnID)
  x_test <- data450

  # reduce x_train and user data (x_test) to common CpG sites
  common <- intersect(colnames(x_train), colnames(x_test))
  x_train <- x_train[,common]
  x_test <- x_test[,common]

  # run models
  if (!is.null(sigma)) { # user supplied a sigma
    rbfker <- kernlab::rbfdot(sigma = sigma)
    y_pred <- metImpute::kernel_pred(x_train,y_train,x_test,ker=rbfker,alpha=0.99)
  }
  else{ # user did not supply a sigma
    sigma <- sigmas[chrom_val, 2] # retrieve sigma associated with chromosome
    rbfker <- kernlab::rbfdot(sigma = sigma)
    y_pred <- metImpute::kernel_pred(x_train,y_train,x_test,ker=rbfker,alpha=0.99)
  }

  # save and export imputed EPIC probes to specified file location
  EPIC <- y_pred
  EPIC_save_path <- paste(path,'/chr',chrom,'_EPIC.rds',sep='')
  saveRDS(EPIC, EPIC_save_path)

  # save and export user data (x_test) with missing values imputed to specified file location
  data450_save_path <- paste(path,'/chr',chrom,'_data450.rds',sep='')
  saveRDS(x_test, data450_save_path)

  # save and export summary statistics table
  summary_stats_save_path <- paste(path, '/chr', chrom, '_summarystats.rds',sep='')
  saveRDS(metImpute::summary_stats, summary_stats_save_path)
}
