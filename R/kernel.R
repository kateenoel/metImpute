#' Kernel
#'
#' Description of function
#' @param data A 450k feature data set with CpG sites as columns
#' @return
#' @export
kernel <- function(data,chrom) {
  # load in training data
  x_train <- get(noquote(paste('chr',chrom,'_xtrain', sep = "")), envir = asNamespace('missingmethyl'), inherits = FALSE)
  y_train <- get(noquote(paste('chr',chrom,'_ytrain', sep = "")), envir = asNamespace('missingmethyl'), inherits = FALSE)

  # impute missing values (source data and user data)


  # reduce x_train and user data to common CpG sites
  common <- intersect(colnames(x_train), colnames(data))
  data <- data[,common]
  x_train <- x_train[,common]

  # run models
  rbfker <- kernlab::rbfdot(sigma = sig)
  y_pred <- kernel_pred(x_train,y_train,x_test,ker=rbfker,alpha=0.99)

  # return 450k feature set + imputed EPIC probes

}
