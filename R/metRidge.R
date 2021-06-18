#' metRidge
#'
#' Use penalized ridge regression to impute level of methylation for EPIC platform CpG sites from an HM450 probe set
#' @param data450 An HM450 probe data set with CpG sites as columns
#' @param data850 An EPIC probe data set with CpG sites as columns
#' @param chrom Indicator of chromosome (e.g. 22 or 'X')
#' @param path Output directory for imputed EPIC probe set and summary statistics
#' @examples
#' metRidge(HM450, EPIC, 22, '/home/project/results')
#' metRidge(HM450, EPIC, 'X', '/home/project/results')
#' @export
metRidge <- function(data450, data850, chrom, path) {
  # create training data objects based on chrom input
  lam <-  get(noquote(paste('chr',chrom,'_lambdas', sep = '')), envir = asNamespace('metImpute'), inherits = FALSE)

  # impute missing values with the column mean (source data and user data)
  for(i in 1:ncol(data850)) {data850[is.na(data850[,i]), i] <- mean(data850[,i], na.rm = TRUE)}

  for(i in 1:ncol(data450)) {data450[is.na(data450[,i]), i] <- mean(data450[,i], na.rm = TRUE)}

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

  # run imputation models
  y_pred <- matrix(NA,nrow=nrow(data),ncol=ncol(y_train)) # set up results matrix
  colnames(y_pred) <- colnames(y_train) # preserve column names
  rownames(y_pred) <- rownames(y_train) # preserve row names
  for(i in 1:ncol(y_train)){
    if(!any(is.na(y_train[,i]))){
      fit <- glmnet::glmnet(x_train, y_train[,i], alpha=0, lambda = lam[i,2])
      y_pred[,i] <- predict(fit, newx = x_test)
    }
  }

  # save and export imputed EPIC probes to specified file location
  EPIC <- y_pred
  EPIC_save_path <- paste(path,'/chr',chrom,'_EPIC.rds',sep='')
  saveRDS(EPIC, EPIC_save_path)

  # save and export summary statistics table
  summary_stats_save_path <- paste(path, '/chr', chrom, '_summarystats.rds',sep='')
  saveRDS(metImpute::summary_stats, summary_stats_save_path)
}
