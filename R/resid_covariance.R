#' @title Residuals Covariance Matrix
#'
#' @description Given an array of matrices, it calculates the residual covariance matrix between column i=1,...,K
#' of the first matrix of the array against column(s) $i=1,...,K of 2,..,k>=2$ matrices of the same array.
#' @param data_array the array containing the different matrices
#' @param k the last column considered for calculating residuals
#' @return the residual covariance matrix
#' @importFrom stats lm
resid_covariance<-function(data_array,k){
  K=ncol(data_array)
  if(k<2){
    stop("k needs to be >=2")
  }
  else{
    res<-sapply(seq_len(K), function(i){lm(data_array[,i,1]~data_array[,i,2:k])$residuals})
  }
  Covar_matr<-(t(res)%*%(res))/(nrow(data_array))
  return(Covar_matr)
}
