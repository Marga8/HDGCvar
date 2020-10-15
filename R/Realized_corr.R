#' @title Computing Realized Correlation
#' @description This function computes realized correlations from realized variances and covariances with the possibility of
#' Fisher-transforming the realized correlations.
#' @param realized_variances Dataset of realized volatilities. A matrix or something that can be coerced to a matrix. Note: the volatilities must not be in logs.
#' @param realized_covariances Dataset of realized covariances. A matrix or something that can be coerced to a matrix.
#' @param fisher_transf Logical: if TRUE the correlations are computed and Fisher-transformed
#' @return LM test statistics and p-values: asymptotic, with finite sample correction and asymptotic with heteroscedasticity correction and Lasso selections are printed to the console
#' @export
#' @examples \dontrun{Realized_corr(real_var, real_cov, fisher_transf=T)}
Realized_corr<-function(realized_variances, realized_covariances,fisher_transf=T){
  Rcov10 = as.matrix(realized_covariances)
  var10= as.matrix(realized_variances)
  if(ncol(Rcov10)!=(((ncol(realized_variances)^2)-ncol(realized_variances))/2)){
    stop(paste("The number of covariances in realized_covariances should be", (((ncol(realized_variances)^2)-ncol(realized_variances))/2),sep=" ") )
  }
  if(log==F){
    var10=exp(var10)#undo the log
    realized_variances=exp(realized_variances)
  }
  stack_v = matrix(NA,nrow(Rcov10),ncol(Rcov10))
  #### Compute correlations ####
  for (i in 1:nrow(var10)) {
    var1<-matrix(var10[i,])
    var2<-c(var10[i,])
    varianze<-as.matrix((var1%*%(var2)))
    varianze[upper.tri(varianze,diag=TRUE)] <-0
    diag(varianze)=1
    dim(varianze)<-c((ncol(var10)*ncol(var10)),1)
    varianze<- varianze[varianze != "0"]
    varianze<- varianze[varianze != "1"]
    stack_v[i,]<-sqrt(varianze)
  }
  Correlations<-Rcov10/stack_v
  ### Fisher-transform the correlations ###
  if(fisher_transf==T){
    realized_correlations = 0.5*(log((1+Correlations)/(1-Correlations)))
  }
  if(fisher_transf==F){
    realized_correlations<-Correlations
  }
  return(realized_correlations)
}
