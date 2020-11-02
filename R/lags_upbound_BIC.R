#' @title Lag length Selection via BIC empirical upper bound
#' @description Selects the lag length p of the VAR using an empirical upper bound: residuals of
#' the diagonalized VAR are used to build the empirical covariance matrix and an approximation of
#' its determinant that uses the matrix trace is employed to be able to select p using Bayesian Information Criterion
#' @param data a dataframe or matrix of the original set of time series forming the VAR
#' @param p_max maximum lag length to consider, default is 10
#' @return  returns the estimated lag length upper bound
#' @export
#' @examples  lags_upbound_BIC(sample_dataset_I1, p_max=10)
#' @references Hecq, A., Margaritella, L., Smeekes, S., "Inference in Non Stationary High Dimensional VARs" (2020, check the latest version at https://sites.google.com/view/luca-margaritella )
#' @references Hecq, A., Margaritella, L., Smeekes, S., "Granger Causality Testing in High-Dimensional VARs: a Post-Double-Selection Procedure." arXiv preprint arXiv:1902.10991 (2019).
lags_upbound_BIC<-function(data,p_max=10){

  data<-as.matrix(data) #data
  K <- ncol(data) #numb of variables
  datalags<-create_lags(data, p = p_max, include.original = TRUE, trim = TRUE) #create p_max lags
  mat_arrays<-split_matrix(datalags,nrow(datalags),K) #split datalags in p_max+1 array's elements

  MatriX_Decision<-matrix(NA,nrow=p_max,ncol=2)
  colnames(MatriX_Decision)<-c("p","BIC")
  MatriX_Decision[,1]<-1:p_max

  for (i in 2:(p_max+1)) {
    Omega<-resid_covariance(mat_arrays,i)
    BIC<-log(prod(diag(Omega)))+((log(nrow(data))/(nrow(data)))*(ncol(data)*(i-1)))
    MatriX_Decision[i-1,2]<-BIC
  }
  MatriX_Decision<-as.data.frame(MatriX_Decision)
  Best<-MatriX_Decision[which.min(MatriX_Decision$BIC),]
  #print(paste("estimated lag-length with BIC is p=",BBB[[1]],sep=""))
  return(Best[[1]])
}
