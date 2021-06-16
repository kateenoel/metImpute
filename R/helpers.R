train_psi_mat <- function(x, ker){
  kermat <- kernlab::kernelMatrix(ker, x)
  n <- ncol(kermat)
  eigres <- eigen(kermat)
  psi_mat <- eigres$vectors
  return(list(psi_mat=psi_mat,lambda=eigres$values,ker=ker))
}

pred_psi_mat <- function(xp,x_train,psi_mat,lambda,ker){
  kp <- kernlab::kernelMatrix(ker,xp,x_train)
  psi_p <- kp%*%psi_mat/matrix(lambda,nrow=nrow(kp),ncol=length(lambda),byrow=TRUE)
  return(psi_p)
}

scatter_pred <- function(y_pred,y_test,idx = 1,pch=19,cex=0.5,main=NULL,ylim=NULL,xlim=NULL){
  plot(y_pred[,1],y_test[,1],ylim=ylim,xlim=xlim,
       pch=pch,xlab="Prediction",ylab="Observation",cex=cex,main=main)
  if(length(idx)>1)
    for(i in idx){
      points(y_pred[,i],y_test[,i],pch=pch,cex=cex,col=i)
    }
  abline(a=0,b=1)
}

rbf_kernel_pred <- function(x_train,y_train,x_test,sigma1 = 0.0001){
  rbf <- rbfdot(sigma = sigma1)
  train_fit <- metImpute::train_psi_mat(x_train,rbf)
  test_psi_mat <- metImpute::pred_psi_mat(x_test,x_train,train_fit$psi_mat,lambda=train_fit$lambda,ker=train_fit$ker)
  beta_coef <- t(train_fit$psi_mat)%*%y_train
  y_pred <- test_psi_mat%*%beta_coef
  return(y_pred)
}

kernel_pred <- function(x_train,y_train,x_test,ker,alpha=0.95,output=NULL){
  train_fit <- metImpute::train_psi_mat(x_train,ker)
  contr <- cumsum(train_fit$lambda)/sum(train_fit$lambda)
  idx <- min(which(contr>alpha))
  test_psi_mat <- metImpute::pred_psi_mat(x_test,x_train,train_fit$psi_mat[,1:idx],
                              lambda=train_fit$lambda[1:idx],ker=train_fit$ker)
  beta_coef <- t(train_fit$psi_mat[,1:idx])%*%y_train
  y_pred <- test_psi_mat%*%beta_coef
  if(is.null(output)){
    return(y_pred)
  }
  return(list(y_pred=y_pred,beta_coef=beta_coef,train_psi_mat=train_fit$psi_mat,
              test_psi_mat=test_psi_mat,lambda=train_fit$lambda,contr=contr,contr_idx=idx))
}

summarize_res <- function(y_pred,y_test){
  rmse <- metImpute::sapply_pb(1:ncol(y_test),function(i) sqrt(mean((y_pred[,i]-y_test[,i])^2,na.rm=TRUE)))
  R <- metImpute::sapply_pb(1:ncol(y_test),function(i) cor(y_pred[,i],y_test[,i],use="pairwise.complete.obs"))
  R2 <- R*R
  mae <- metImpute::sapply_pb(1:ncol(y_test), function(i) mean(abs(y_pred[,i]-y_test[,i])))
  tab <- c(rmse=mean(rmse,na.rm=TRUE),R2=mean(R2,na.rm=TRUE), mae=mean(mae, na.rm=TRUE))
  return(list(rmse=rmse,R2=R2,mae=mae,tab=tab))
}

sapply_pb <- function(X, FUN, ...){
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- sapply(X, wrapper, ...)
  close(pb)
  res
}

