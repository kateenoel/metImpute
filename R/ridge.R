#' Ridge regression
#'
#' Description of function
#' @param data A 450k feature data set with CpG sites as columns
#' @param chrom Indicator of chromosome
#' @return
#' @export
ridge <- function(data, chrom) {
  # load required source files (training data, functions, lambdas, subset indicators)


  # impute missing values (source data and user data)


  # reduce x_train and user data to common CpG sites
  common <- intersect(colnames(x_train), colnames(data))
  data <- data[,common]
  x_train <- x_train[,common]

  # run imputation models
  y_pred <- matrix(NA,nrow=nrow(x_test),ncol=ncol(y_train)) # set up results matrix
  colnames(y_pred) <- colnames(y_train) # preserve column names
  for(i in 1:ncol(y_train)){
    if(!any(is.na(y_train[,i]))){
      fit=glmnet(x_train,y_train[,i], alpha=0, lambda = lam[i,2])
      y_pred[,i] <- predict(fit,newx=x_test)
    }
  }

  # return

}
