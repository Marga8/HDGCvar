#' @title Test Granger causality in High Dimensional Stationary VARs
#'
#' @param  GCpair     a named list with names GCto and GCfrom containing vectors of the relevant GC variables.
#' @param  data       data matrix or object that can be coerced to a matrix.
#' @param  p          lag length of the VAR
#' @param  bound      lower bound on tuning parameter lambda
#' @param  parallel   TRUE for parallel computing
#' @param  n_cores    nr of cores to use in parallel computing, default is all but one
#' @return            LM Chi-square test statistics (asymptotic), LM F-stat with finite sample correction, both with their corresponding p-value.
#' Lasso selections are also printed to the console.
#' @export
#' @importFrom parallel makeCluster clusterSetRNGStream clusterExport clusterEvalQ detectCores parSapply stopCluster
#' @importFrom stats cor
#' @examples HDGC_VAR_I0(GCpair=list("GCto"="Var 1", "GCfrom"="Var 2"), data=sample_dataset_I0, p=2)
HDGC_VAR_I0 <- function(GCpair, data, p = 1, bound = 0.5 * nrow(data),
                        parallel = FALSE, n_cores = NULL) {
  GCto <- GCpair$GCto #Granger-caused variable
  GCfrom <- GCpair$GCfrom #Granger-causing
  data <- scale(as.matrix(data), scale = FALSE) #scale dataset
  Mat_corr <- as.matrix(cor(data))
  upper <- as.vector(Mat_corr[upper.tri(Mat_corr)])
  amount <- sum(upper >(1-0.001))
  if(any(1-0.001<upper)){
    warning(paste("The used dataset contains",amount, "correlations larger than 0.999, this can cause failure of OLS and poor variable selection.",sep=" "))
  }
  K <- ncol(data) #numb of variables
  X_all <- create_lags(data, p , include.original = FALSE)  #create p lags + d augmentation of all (original K not included) #correspond to lDatafin=xcont
  if((ncol(X_all)+ncol(data))>nrow(data)){
    warning( paste("You are estimating an HD model in which each equation has p*K=",(ncol(X_all)+ncol(data)), " parameters and ", nrow(data), " observations.
                   Depending on how large is p, to avoid failure of OLS you might want to decrease the bound=0.5*nrow(data)." ))
  }
  X <- X_all[, 1:(K * p)] #original K*p lags of regressors , correspond to lDatafin1
  Y <- data[-(1:p), ] #original variables , correspond to OriginalV
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
  GCfrom_names<-colnames(X)[X_index] #names of the p GCfrom lags
  GCto_names<-colnames(X)[Y_index] #names of the p GCto lags
  X_GC <- X[, X_index] #corresponds to dataframe3
  Z <- X[, -X_index] #all other variables but Granger-causing, corresponds to xcontnogc
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
    clusterEvalQ(cl = cl, library(glmnet))

    lasso_Sx <- parSapply(cl, seq_len(length(X_index)),
                          active_set_1, d=0,p=p, y = X_GC, X_index =NULL, z = Z, z_a =NULL , bound = bound) #regress each time eliminating only the corresponding granger causing lag
    stopCluster(cl)
  }
  else {
    lasso_Sx <- sapply(seq_len(length(X_index)),
                       active_set_1, d=0,p=p, y = X_GC, X_index =NULL, z = Z, z_a =NULL , bound = bound) #regress each time eliminating only the corresponding granger causing lag
  }

  # Regression for y
  lasso_Sy <- active_set_1(i = 1, d = 0, p = p, X_index = NULL, y = y_I,
                           z = Zx, z_a = NULL, bound = bound)

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
  LM_out <- LM_test(y_I, diag(I) %x% X_GC, Zx_sel, I) #dep variable, Granger causing original p lags, selected variables
  return(list(tests = LM_out, selections = lasso_sel))
}
