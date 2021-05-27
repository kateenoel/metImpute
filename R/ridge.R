#' Ridge regression
#'
#' Description of function
#' @param data A 450k feature data set with CpG sites as columns
#' @param chrom Indicator of chromosome
#' @return
#' @export
ridge <- function(data, chrom) {
  # create training data objects based on chrom input
  x_train <- get(noquote(paste('chr',chrom,'_xtrain', sep = "")), envir = asNamespace('missingmethyl'), inherits = FALSE)
  y_train <- get(noquote(paste('chr',chrom,'_ytrain', sep = "")), envir = asNamespace('missingmethyl'), inherits = FALSE)
  lam <-  get(noquote(paste('chr',chrom,'_lambdas', sep = "")), envir = asNamespace('missingmethyl'), inherits = FALSE)

  # impute missing values (source data and user data)


  # reduce x_train and user data to common CpG sites
  common <- intersect(colnames(x_train), colnames(data))
  data <- data[,common]
  x_train <- x_train[,common]

  # run imputation models
  y_pred <- matrix(NA,nrow=nrow(data),ncol=ncol(y_train)) # set up results matrix
  colnames(y_pred) <- colnames(y_train) # preserve column names
  for(i in 1:ncol(y_train)){
    if(!any(is.na(y_train[,i]))){
      fit=glmnet(x_train,y_train[,i], alpha=0, lambda = lam[i,2])
      y_pred[,i] <- predict(fit,newx=data)
    }
  }

  # return

}
