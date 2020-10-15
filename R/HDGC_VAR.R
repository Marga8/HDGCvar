#' @title  Test Granger causality in High Dimensional mixed Integrated and Cointegrated VARs
#'
#' @param  GCpair     a named list with names GCto and GCfrom containing vectors of the relevant GC variables.
#' @param  data       a data matrix or something that can be coerced to a matrix
#' @param  p          lag length of the VAR
#' @param  d          order of lag augmentation corresponding to suspected max order of integration
#' @param  bound      lower bound on tuning parameter lambda
#' @param  parallel   TRUE for parallel computing
#' @param  n_cores    nr of cores to use in parallel computing, default is all but one
#' @return            LM test statistics, p-values: asymptotic and with finite sample correction and Lasso selections are printed to the console
#' @export
#' @importFrom parallel makeCluster clusterSetRNGStream clusterExport clusterEvalQ detectCores parSapply stopCluster parLapply
#' @importFrom stats cor
#' @examples \dontrun{GCpair<-list("GCto"="X1", "GCfrom"="X2")
#' HDGC_VAR(GCpair, data, p=2,d=2, parallel=T)}
HDGC_VAR <- function(GCpair, data, p = 1, d = 0, bound = 0.5 * nrow(data),
                     parallel = FALSE, n_cores = NULL) {
  GCto <- GCpair$GCto #Granger-caused variable
  GCfrom <- GCpair$GCfrom #Granger-causing
  data<-as.matrix(data)
  Mat_corr<-as.matrix(cor(data))
  upper<-as.vector(Mat_corr[upper.tri(Mat_corr)])
  amount<-sum(upper >(1-0.001))
  if(any(1-0.001<upper)){
    warning(paste("The used dataset contains",amount, "correlations larger than 0.999, this can cause failure of OLS and poor variable selection.",sep=" "))
  }
  K <- ncol(data) #numb of variables
  X_all <- create_lags(data, p + d, include.original = FALSE)  #create p lags + d augmentation of all (original K not included) #correspond to lDatafin=xcont
  if((ncol(X_all)+ncol(data))>nrow(data)){
    warning( paste("You are estimating an HD model in which each equation has p*K=",(ncol(X_all)+ncol(data)), " parameters and ", ncol(data), " observations.
                   Depending on how large is p, to avoid failure of OLS you might want to increase the bound=0.5*nrow(data)." ))
  }
  X <- as.matrix(X_all[, 1:(K * p)]) #original K*p lags of regressors , correspond to lDatafin1
  Y <- as.matrix(data[-(1:(p + d)), ]) #original variables , correspond to OriginalV, trimmed of the lags NA
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
  ## added
  Y_index<-c(sapply(y_index, seq, by = K, length.out = p, simplify = "array"))
  ## added
  GCfrom_names<-colnames(X)[X_index] #names of the p GCfrom lags
  GCto_names<-colnames(X)[Y_index]

  X_GC <- as.matrix(X[, X_index]) #corresponds to dataframe3
  Z <- X[, -X_index] #all other variables but Granger-causing, corresponds to xcontnogc

  ##added
  Z_index<-match( GCto_names, colnames(Z) ) #indeces of GCto_names after removinf p lags of GCfrom


  ZX_augm <- NULL
  if (d > 0) {
    X_augm <- X_all[, K * p + 1:(K * d)] #lags p+1 ..p+d of all variables, i.e all extra lags
    Z_augm <- X_augm[, c(sapply(c(y_index, x_index), seq, by = K, length.out = d,
                                simplify = "array"))]# lags p+1 ..p+d of Granger caused and Granger causing
    #if (d + 1 > p) {
    # ZX_augm <- X_augm[, c(sapply(c(y_index, x_index), seq, by = K, length.out = d + 1 - p,
    #                             simplify = "array"))]
    #}
  }
  Zx <- diag(I) %x% X #full matrix including p lags of the grangercausing
  Zx_1 <- diag(I) %x% Z ## matrix without lags of the grangercausing
  # Regressions for X_GC
  if (parallel) {
    if (is.null(n_cores)) {
      n_cores <- detectCores() - 1
    }
    cl <- makeCluster(n_cores, setup_strategy = "sequential")
    clusterSetRNGStream(cl, sample.int(2^20, size = 1))
    clusterExport(cl = cl, Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv)))
    clusterEvalQ(cl = cl, library(glmnet))

    lasso_Sx_1 <- parLapply(cl, seq_len(length(X_index)),
                            active_set_1, d=d, p=p, y = X_GC, X_index =X_index, z = X, z_a =NULL , bound = bound) #regress each time eliminating only the corresponding granger causing lag
    stopCluster(cl)
  }
  else {
    lasso_Sx_1 <- lapply(seq_len(length(X_index)),
                         active_set_1, d=d, p=p, y = X_GC, X_index =X_index, z = X, z_a =NULL , bound = bound) #regress each time eliminating only the corresponding granger causing lag
  }
  if (p==1){ #when p=1 there's no need to exclude the Granger causing lags as active_set already exclude the single relevant lag
    lasso_Sx<-as.matrix(lasso_Sx_1[[1]])
  }
  if (p>1){
    lasso_Sx<-sapply(seq_len(length(X_index)), function(i){
      as.matrix(lasso_Sx_1[[i]])[-which(names(lasso_Sx_1[[i]]) %in% GCfrom_names),] })
  }
  # Regression for y
  lasso_Sy_1 <- active_set_1(i = 1,d=d, p=p, X_index=NULL , y = y_I, z = Zx, z_a = NULL, bound = bound)
  lasso_Sy<-as.matrix(lasso_Sy_1[-c(X_index)]) # cut off the elements corresponding to GC lags
  rownames(lasso_Sy)<-c(1:nrow(lasso_Sy))

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
    Zx_sel <- Zx_1[, lasso_sel]
  }

  # Construct lag augmentations
  if (d > 0) { #augmentation of x and y
    Zx_sel <- cbind(Zx_sel, diag(I) %x% Z_augm)
  }

  # Perform LM test
  LM_out <- LM_test(y_I, diag(I) %x% X_GC, Zx_sel, I) #dep variable, Granger causing original p lags, selected variables
  return(list(tests = LM_out, selections = lasso_sel))
}
