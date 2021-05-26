#' Kernel
#'
#' Description of function
#' @param data A 450k feature data set with CpG sites as columns
#' @return
#' @export
kernel <- function(data) {
  # load required source files (training data, functions, subset indicators)


  # impute missing values (source data and user data)


  # subset training data into 450k feature set and EPIC-only feature set


  # reduce x_train and user data to common CpG sites
  common <- intersect(colnames(x_train), colnames(data))
  data <- data[,common]
  x_train <- x_train[,common]

  # run models
  rbfker <- kernlab::rbfdot(sigma = sig)
  y_pred <- kernel_pred(x_train,y_train,x_test,ker=rbfker,alpha=0.99)

  # return 450k feature set + imputed EPIC probes

}
