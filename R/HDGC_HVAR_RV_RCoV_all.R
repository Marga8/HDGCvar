#' @title Networks of Realized Volatilities conditional on the set of Realized Correlations
#'
#' @param  realized_variances        Dataset of (stationary) realized volatilities. A matrix or object that can be coerced to a matrix.
#' @param  realized_covariances      Dataset of (stationary) realized covariances. A matrix or object that can be coerced to a matrix. Note: the columns should exactly
#' be (((ncol(realized_volatilities)^2)-ncol(realized_volatilities))/2)
#' @param  fisher_transf             Logical: if TRUE the correlations are computed and Fisher transformed
#' @param  log                       Default is TRUE, if the realized volatilities are already log transformed then put to FALSE
#' @param  bound                     Lower bound on lambda
#' @param  parallel                  TRUE for parallel computing
#' @param  n_cores                   Nr of cores to use in parallel computing, default is all but one
#' @return   Granger causality matrix and Lasso selections are printed to the console
#' @export
#' @examples \dontrun{ HDGC_HVAR_RV_RCoV_all(real_var, real_cov, fisher_transf=T, log=TRUE ,parallel = TRUE) }
#' @references Hecq, A., Margaritella, L., Smeekes, S., "Granger Causality Testing in High-Dimensional VARs: a Post-Double-Selection Procedure." arXiv preprint arXiv:1902.10991 (2019).
#' @references  Corsi, Fulvio. "A simple approximate long-memory model of realized volatility." Journal of Financial Econometrics 7.2 (2009): 174-196.
HDGC_HVAR_RV_RCoV_all <- function(realized_variances,realized_covariances, fisher_transf=TRUE, log=TRUE, bound = 0.5 * nrow(realized_variances),
                                  parallel = FALSE, n_cores = NULL) {

  #### Load realized variances and covariances ####
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
  # take log of realized variances
  realized_variances<-log(realized_variances)

  varnames <- colnames(realized_variances)
  K <- ncol(realized_variances)
  GCpairs <- vector("list", length = K * (K - 1))
  ind <- 0
  for (i in 1:K) {
    for (j in (1:K)[-i]) {
      ind <- ind + 1
      GCpairs[[ind]] <- list(GCto = varnames[i], GCfrom = varnames[j])
    }
  }
  GC_all_pairs <- HDGC_HVAR_multiple_RVCOV(realized_variances=realized_variances, realized_correlations=realized_correlations,
                                           GCpairs = GCpairs, log=log, bound = bound, parallel = parallel, n_cores = n_cores)
  GC_matrix <- array(dim = c(K, K, 2, 3))
  dimnames(GC_matrix) <- list(GCto = varnames, GCfrom = varnames,
                              stat = c("LM_stat", "p_value"), type = c("Asymp", "FS_cor", "Asymp_Robust"))
  for (i in 1:length(GCpairs)) {
    ind_to <- which(varnames %in% GCpairs[[i]]$GCto)
    ind_from <- which(varnames %in% GCpairs[[i]]$GCfrom)
    GC_matrix[ind_to, ind_from, , ] <- GC_all_pairs$tests[, , i]
  }
  return(list(tests = GC_matrix, selections = GC_all_pairs$selections))
}
