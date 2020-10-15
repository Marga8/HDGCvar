#' @title  Test Granger causality for Realized Volatilities in High Dimensional Stationary Heterogeneous VARs conditioning on Realized Correlations
#'
#' @param  GCpair     A named list with names GCto and GCfrom containing vectors of the relevant GC variables.
#' @param realized_variances Dataset of realized volatilities. A matrix or something that can be coerced to a matrix. Note: the volatilities must not be in logs.
#' @param realized_correlations Dataset of realized correlations. To compute realized correlations from realized variances and realized covariances use \code{\link{Realized_corr}}
#' @param  bound      lower bound on tuning parameter lambda
#' @param  parallel   TRUE for parallel computing
#' @param  n_cores    nr of cores to use in parallel computing, default is all but one
#' @return            LM test statistics and p-values: asymptotic, with finite sample correction and asymptotic with heteroscedasticity correction and Lasso selections are printed to the console
#' @export
#' @importFrom parallel makeCluster clusterSetRNGStream clusterExport clusterEvalQ detectCores parSapply stopCluster parLapply
#' @examples \dontrun{GCpair<-list("GCto"="X", "GCfrom"="Z")
#' HDGC_HVAR_RVCOV(GCpair, real_var, real_corr,parallel = T)}
HDGC_HVAR_RVCOV <- function(GCpair, realized_variances, realized_correlations, bound = 0.5 * nrow(realized_variances),
                            parallel = FALSE, n_cores = NULL) {
  p = 3 #impose Three lags (HVAR)

  GCto <- GCpair$GCto #Granger-caused variable
  GCfrom <- GCpair$GCfrom #Granger-causing
  realized_variances<-as.matrix(realized_variances)
  realized_correlations<-as.matrix(realized_correlations)

  K_1<-ncol(realized_variances)
  K_2<-ncol(realized_correlations)
  K <- K_1+K_2 # tot numb of variables

  databind<-cbind(realized_variances,realized_correlations)
  X_all<-create_lags_RV(databind,include.original = FALSE) #all lags
  X_all<-as.matrix(X_all)

  Y <- realized_variances[-(1:22), ] #original variables cut the same burn in
  realized_correlations<-realized_correlations[-(1:22),] #original variables cut the same burn in

  y_index <- which(colnames(Y) %in% GCto) #index of Granger-caused variable
  if (is.null(y_index)) {
    stop("No matching variable for GCto found.")
  }
  I <- length(y_index) #number of dep variables
  y_I <- c(Y[, y_index]) #dependent variable, corresponds to ycont1
  x_index <- which(colnames(Y) %in% GCfrom) #index of Granger-causing variable
  if (is.null(x_index)) {
    stop("No matching variable for GCfrom found.")
  }
  X_index <- c(sapply(x_index, seq, by = K, length.out = p, simplify = "array")) #indeces of p lags of Granger causing
  Y_index<-c(sapply(y_index, seq, by = K, length.out = p, simplify = "array"))
  GCfrom_names<-colnames(X_all)[X_index] #names of the p GCfrom lags
  GCto_names<-colnames(X_all)[Y_index] #names of the p GCto lags

  X_GC <- X_all[, X_index] #corresponds to dataframe3
  Z <- X_all[, -X_index] #all other variables but Granger-causing, corresponds to xcontnogc

  Z_index<-match( GCto_names, colnames(Z) ) #indeces of GCto_names after removing p lags of GCfrom

  Zx <- diag(I) %x% Z   ## takes off the lags of the grangercausing, corresponds to xcontnogc

  # Regressions for X_GC
  if (parallel) {
    if (is.null(n_cores)) {
      n_cores <- detectCores() - 1
    }
    cl <- makeCluster(n_cores, setup_strategy = "sequential")
    clusterSetRNGStream(cl, sample.int(2^20, size = 1))
    clusterExport(cl = cl, Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv)))
    clusterEvalQ(cl = cl, {
      library(glmnet)
      library(zoo)
    })

    lasso_Sx <- parSapply(cl, seq_len(length(X_index)),
                          active_set_1, d=0,p=p, y = X_GC, X_index =NULL, z = Z, z_a =NULL , bound = 0.5*nrow(Z))
    stopCluster(cl)
  }
  else {
    lasso_Sx <- sapply(seq_len(length(X_index)),
                       active_set_1, d=0,p=p, y = X_GC, X_index =NULL, z = Z, z_a =NULL , bound = 0.5*nrow(Z))
  }

  # Regression for y
  lasso_Sy <- active_set_1(i = 1, d=0,p=p, X_index=NULL , y = y_I, z = Z, z_a = NULL, bound = 0.5*nrow(Z))

  # Collect all active sets
  lasso_S <- cbind((rep(1, I) %x% lasso_Sx) == TRUE, lasso_Sy)

  # Force lags of dependent variable inside the union
  lasso_S[Z_index,]<-c(rep(T,length(X_index)+1)) #in case lasso has kicked out (i.e. put to false) lags of GCto, we should force them in; +1 is for the constant

  # Union of selected variables
  lasso_sel <- apply(lasso_S, 1, any)

  if (I == 1) {
    names(lasso_sel) <- colnames(Z)
  } else {
    names(lasso_sel) <- outer(1:I, colnames(Z), paste, sep = "_")
  }
  if (!any(lasso_sel)) {
    Zx_sel <- NULL
  } else {
    Zx_sel <- Zx[, lasso_sel]
  }

  # Perform LM test
  LM_out <- LM_test_robust(y_I, diag(I) %x% X_GC, Zx_sel, I) #dep variable, Granger causing original p lags, selected variables
  return(list(tests = LM_out, selections = lasso_sel))
}
