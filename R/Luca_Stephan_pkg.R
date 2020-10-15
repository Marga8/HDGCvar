#######################################################################################################
####################################### AUXILIARY FUNCTIONS ###########################################
#' @title Lags creation
#'
#' @description Creates a lag matrix of order p. Can include or exclude original series, and trim the NAs in the
#' start of the sample.
#' @param y vector or matrix to lag
#' @param p lag length
#' @param include.original logical, if TRUE the original series are left inside the output matrix.
#' @param trim logical, if TRUE the initial NAs due to the lag gets trimmed
#' @return the lagged vector or matrix
create_lags <- function(y, p = 1, include.original = TRUE, trim = TRUE) {
  x <- as.matrix(y)
  n <- nrow(x)
  k <- ncol(x)
  lx <- matrix(0, nrow = n, ncol = (p + include.original) * k)
  if (is.null(colnames(x))) {
    c.names <- rep("", k)
  } else {
    c.names <- colnames(x)
  }
  colnames(lx) <- rep(c.names, p + include.original)
  for (i in (1 - include.original):p) {
    cols <- k * (i - 1 + include.original) + 1:k
    lx[(1+i):n, cols] <- x[1:(n-i), ]
    colnames(lx)[cols] <- paste(c.names, " l", i, sep = "")
  }
  return(lx[(1 + trim * p):n, , drop = FALSE])
}

#' @title Lags creation: Daily, Weekly, Monthly aggregation for Realized Volatilities
#'
#' @description Creates a matrix of order 3 containing Daily, Weekly, Monthly returns for Realized Volatilities.
#' Can include or exclude original series, and trim the NAs in the start of the sample. It should be used for HVARs.
#' @param y vector or matrix to lag
#' @param include.original logical, if TRUE the original series are left inside the output matrix.
#' @param trim logical, if TRUE the initial NAs due to the lag gets trimmed
#' @return the lagged vector or matrix
#' @importFrom zoo rollapply
create_lags_RV <- function(y, include.original = TRUE, trim = TRUE) {
  x <- as.matrix(y)
  n <- nrow(x)
  k <- ncol(x)
  lx <- matrix(0, nrow = n, ncol = 4 * k)
  if (is.null(colnames(x))) {
    c.names <- rep("", k)
  } else {
    c.names <- colnames(x)
  }
  colnames(lx) <- rep(c.names, 4)

  for (i in 0:0) {
    cols <-  1:k
    lx[(1+i):n, cols] <- x[1:(n), ]
    colnames(lx)[cols] <- colnames(x)
  }
  for (i in 1:1) {
    cols <- (k+1):(2*k)
    lx[(1+i):n, cols] <- x[1:(n-i), ]
    colnames(lx)[cols] <- paste(c.names, " Daily", " lag", sep = "")
  }
  for (i in 2:2) {
    cols <- ((2*k)+1):(3*k)
    lx[1:5,cols]<-0
    lx[(6):n, cols] <- c(rollapply(x[1:(n-1), ], width = 5, by = 1, FUN = mean, align = "left"))
    colnames(lx)[cols] <- paste(c.names, " Weekly", " lag", sep = "")
  }
  for (i in 3:3) {
    cols <- ((3*k)+1):(4*k)
    lx[1:22, cols]<-0
    lx[(23):n, cols] <- rollapply(x[1:(n-1), ], width = 22, by = 1, FUN = mean, align = "left")
    colnames(lx)[cols] <- paste(c.names, " Monthly", " lag", sep = "")
  }
  if(include.original==F){
    lx<-lx[,(k+1):(4*k)]
  }
  return(lx[(1 + trim * 22):n, , drop = FALSE])
}


#' @title Split Big Matrix
#'
#' @description Splits a big matrix in sub matrices using Kronecker product
#' @param M the matrix
#' @param r the rows of desired cut.
#' @param c the columns of desired cut.
#' @return the splitted matrix
split_matrix <- function(M, r, c){
  nr <- ceiling(nrow(M)/r)
  nc <- ceiling(ncol(M)/c)
  newM <- matrix(NA, nr*r, nc*c)
  newM[1:nrow(M), 1:ncol(M)] <- M

  div_k <- kronecker(matrix(seq_len(nr*nc), nr, byrow = TRUE), matrix(1, r, c))
  matlist <- split(newM, div_k)
  N <- length(matlist)
  mats <- unlist(matlist)
  dim(mats)<-c(r, c, N)
  return(mats)
}

#'@title Ordinary Least Squares
#'
#'@description Perform least squares of y on x. Returns a list with the coefficients, t statistics, p values and
#' residual vector.
#' @param y vector dependent variable
#' @param x a vector or matrix of regressors
#' @return coefficients, fitted values, t statistics, p values, residuals
#' @importFrom stats pnorm
ols <- function(y, x) { ## There is NO INTERCEPT
  y <- as.matrix(y)
  #x<-as.matrix(x)
  x <- as.matrix(cbind(int=1,x)) ## added
  xx.inv <- chol2inv(chol(crossprod(x)))
  b <- xx.inv %*% crossprod(x, y)
  fit <- x %*% b
  e <- y - fit
  s2 <- mean(e^2)
  se.mat <- s2 * xx.inv
  t.stats <- b / sqrt(diag(se.mat))
  p.val <- 2*pnorm(abs(t.stats), lower.tail = FALSE)
  return(list(coef = b, fit = fit, tstat = t.stats, pvalue = p.val, resid = e))
}

#'@title Ordinary Least Squares without intercept
#'
#'@description Perform least squares of y on x. Returns a list with the coefficients, t statistics, p values and
#' residual vector.
#' @param y vector dependent variable
#' @param x a vector or matrix of regressors
#' @return coefficients, fitted values, t statistics, p values, residuals
#' @importFrom stats pnorm
ols_no_int<- function(y, x) { ## There is NO INTERCEPT
  y <- as.matrix(y)
  x<-as.matrix(x)
  #x <- as.matrix(cbind(int=1,x)) ## added
  xx.inv <- chol2inv(chol(crossprod(x)))
  b <- xx.inv %*% crossprod(x, y)
  fit <- x %*% b
  e <- y - fit
  s2 <- mean(e^2)
  se.mat <- s2 * xx.inv
  t.stats <- b / sqrt(diag(se.mat))
  p.val <- 2*pnorm(abs(t.stats), lower.tail = FALSE)
  return(list(coef = b, fit = fit, tstat = t.stats, pvalue = p.val, resid = e))
}

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


#' @title Lasso Bayesian Information Criterion
#'
#' @description Do BIC on lasso output. Needs degrees of freedom coming from glmnet. Bound on lambda is optional,
#' imposed as 0.5T by default.
#' @param y vector of dependent variable
#' @param x vector or matrix of regressors
#' @param b lasso coefficients (taken from glmnet)
#' @param df lasso degrees of freedom (taken from glmnet)
#' @param bound upper bound on tuning parameter lambda: i.e. impose lasso to stop the selection after it selects 0.5T parameters.
#' @return the estimated parameters with BIC
bic <- function(y, x, b, df, bound = 0.5 * length(y)) {
  k <- sum(df <= bound)
  n <- length(y)
  bic <- rep(NA,k)
  for (i in 1:k) {
    e <- y - cbind(1,x) %*% b[,i]
    s.2 <- crossprod(e) / n
    bic[i] <- log(s.2) + df[i] * log(n) / n
  }
  b.bic <- b[, which.min(bic)]
  return(b.bic)
}

#' Calculate square root of a matrix
#'
#' @param M the matrix
#' @return the matrix square rooted
matrix_sqrt <- function(M) {
  eva_ve <- eigen(M)
  S <- eva_ve$vectors %*% diag(sqrt(eva_ve$values))
  return(S)
}

#' Inverse square root matrix
#'
#' @param M the matrix
#' @return the inverse matrix square rooted
inv_matrix_sqrt <- function(M) {
  if (NROW(M) > 1) {
    eva_ve <- eigen(M)
    S <- diag(sqrt(1/eva_ve$values)) %*% t(eva_ve$vectors)
  } else {
    S <- 1 / sqrt(M)
  }
  return(S)
}

#'@title Determine the active set variables via lasso regression but it augments also the lasso selection.
#'
#'@param i determines which column of y is considered
#'@param d the decided lag-augmentation of the interest variables
#'@param p the lag length
#'@param X_index index of the Granger causing variables
#'@param y vector/matrix of dependent variable
#'@param z vector/matrix of regressors
#'@param z_a those variables that should not be penalized
#'@param bound upper bound on tuning parameter lambda
#'@return the active set estimated via lasso using BIC as penalty
#'@importFrom glmnet glmnet
#'@importFrom stats coef
active_set <- function(i,d,p, X_index, y, z, z_a = NULL, bound = 0.5 * NROW(z)) {
  if(is.null(X_index)){
    if (is.null(z_a)) {
      x <- z
      w <- rep(1, ncol(x))
    } else {
      x <- cbind(z, z_a[, seq_len(d) + i - p])
      w <- c(rep(1, ncol(z)), 0)
    }
  }
  else {
    if (is.null(z_a)) {
      x <- z[,-X_index[i]]
      w <- rep(1, ncol(x))
    } else {
      x <- cbind(z[,-X_index[i]], z_a[, seq_len(d) + i - p])
      w <- c(rep(1, ncol(z)-1), 0)
    }
  }
  y <- as.matrix(y)
  lasso_out <- glmnet(x, y[, i], penalty.factor = w, standardize=F)
  lasso_betas <- coef(lasso_out)
  lasso_df <- lasso_out$df
  lasso_bic_beta <- bic(y[, i], x, lasso_betas, lasso_df, bound = bound)
  act_set <- abs(lasso_bic_beta[2:length(lasso_bic_beta)]) > 0
  return(act_set[1:(ncol(z)-!is.null(X_index))])
}

#'@title Determine the active set variables via lasso regression.
#'
#'@param i determines which column of y is considered
#'@param d the decided lag-augmentation of the interest variables
#'@param p the lag length
#'@param X_index index of the Granger causing variables
#'@param y vector/matrix of dependent variable
#'@param z vector/matrix of regressors
#'@param z_a those variables that should not be penalized
#'@param bound upper bound on tuning parameter lambda
#'@return the active set estimated via lasso using BIC as penalty
#'@importFrom glmnet glmnet
#'@importFrom stats coef
active_set_1 <- function(i,d,p, X_index, y, z, z_a = NULL, bound = 0.5 * NROW(z)) {
  if(is.null(X_index)){
    if (is.null(z_a)) {
      x <- z
      w <- rep(1, ncol(x))
    } else {
      x <- cbind(z, z_a[, seq_len(d) + i - p])
      w <- c(rep(1, ncol(z)), 0)
    }
  }
  else {
    if (is.null(z_a)) {
      x <- z[,-X_index[i]]
      w <- rep(1, ncol(x))
    } else {
      x <- cbind(z[,-X_index[i]], z_a[, seq_len(d) + i - p])
      w <- c(rep(1, ncol(z)-1), 0)
    }
  }
  y <- as.matrix(y)
  lasso_out <- glmnet(x, y[, i], penalty.factor = w, standardize=F)
  lasso_betas <- coef(lasso_out)
  lasso_df <- lasso_out$df
  lasso_bic_beta <- bic(y[, i], x, lasso_betas, lasso_df, bound = bound)
  act_set <- abs(lasso_bic_beta[2:length(lasso_bic_beta)]) > 0
  return(act_set)
}

#'@title LM test
#'@description Perform the LM test of significance of x for y given z.
#'@param y dependent variable: Granger-caused
#'@param x Granger causing variable, all its lags
#'@param z conditioning set of variables
#'@param I  the numberr of dependent variables, that is number of equations in the VAR.
#'@return test statistics, both asymptotic chi square and finite sample corrected F and relative p values
#'@importFrom stats pchisq pf
LM_test <- function(y, x, z, I = 1) {
  if (is.null(z)) {
    xi_mat <- matrix(y, ncol = I)
    n <- nrow(xi_mat)
    Sigma_hat <- (crossprod(xi_mat) / n) %x% diag(n)

    y_star <- inv_matrix_sqrt(Sigma_hat) %*% y
    x_star <- inv_matrix_sqrt(Sigma_hat) %*% x

    ols2_out <- ols(y_star, x_star)
    nu_star <- ols2_out$resid
    R_sq <- (crossprod(y_star) - crossprod(nu_star)) / nrow(y_star)
    LM_stat1 <- n * R_sq
    p_value1 <- pchisq(LM_stat1, ncol(x_star), lower.tail = FALSE)
    LM_stat2 <- ((nrow(y_star) - ncol(x_star)) / ncol(x_star)) * R_sq / (1 - R_sq)
    p_value2 <- pf(LM_stat2, ncol(x_star), nrow(y_star) - ncol(x_star), lower.tail = FALSE)
  } else {
    ols_out <- ols(y, z)
    xi <- ols_out$resid
    xi_mat <- matrix(xi, ncol = I)
    n <- nrow(xi_mat)
    Sigma_hat <- (crossprod(xi_mat) / n) %x% diag(n)

    Inv_Sigma<-inv_matrix_sqrt(Sigma_hat)
    y_star <- Inv_Sigma %*% y
    z_star <- Inv_Sigma %*% z
    x_star <- Inv_Sigma %*% x

    ols2_out <- ols(y_star, z_star)
    xi_star <- ols2_out$resid
    ols2_out <- ols(xi_star, cbind(x_star, z_star))
    nu_star <- ols2_out$resid

    LM_stat1 <- crossprod(xi_star) - crossprod(nu_star)
    p_value1 <- pchisq(LM_stat1, ncol(x_star), lower.tail = FALSE)
    LM_stat2 <- ((nrow(y_star) - ncol(x_star) - ncol(z_star)) / ncol(x_star)) *
      LM_stat1 / (nrow(y_star) - LM_stat1)
    p_value2 <- pf(LM_stat2, ncol(x_star), nrow(y_star) - ncol(x_star) - ncol(z_star),
                   lower.tail = FALSE)
  }
  LM_out <- matrix(c(LM_stat1, p_value1, LM_stat2, p_value2), nrow = 2, ncol = 2)
  rownames(LM_out) <- c("LM_stat", "p_value")
  colnames(LM_out) <- c("Asymp", "FS_cor")
  return(LM_out)
}

#'@title LM test heteroscedastic robust
#'@description Perform the LM test of significance of x for y given z and add the robust asymptotic test statistic.
#'@param y dependent variable: Granger-caused
#'@param x Granger causing variable, all its lags
#'@param z conditioning set of variables
#'@param I  the numberr of dependent variables, that is number of equations in the VAR.
#'@return test statistics, both asymptotic chi square, finite sample corrected F, asymptotic chi square robust to heteroscedasticity and relative p values
#'@importFrom stats pchisq pf
LM_test_robust <- function(y, x, z, I = 1) {
  if (is.null(z)) {
    xi_mat <- matrix(y, ncol = I)
    n <- nrow(xi_mat)
    Sigma_hat <- (crossprod(xi_mat) / n) %x% diag(n)

    y_star <- inv_matrix_sqrt(Sigma_hat) %*% y
    x_star <- inv_matrix_sqrt(Sigma_hat) %*% x

    ols2_out <- ols(y_star, x_star)
    nu_star <- ols2_out$resid
    R_sq <- (crossprod(y_star) - crossprod(nu_star)) / nrow(y_star)
    LM_stat1 <- n * R_sq
    p_value1 <- pchisq(LM_stat1, ncol(x_star), lower.tail = FALSE)
    LM_stat2 <- ((nrow(y_star) - ncol(x_star)) / ncol(x_star)) * R_sq / (1 - R_sq)
    p_value2 <- pf(LM_stat2, ncol(x_star), nrow(y_star) - ncol(x_star), lower.tail = FALSE)
    LM_stat3<-LM_stat1
    p_value3<-p_value1
  } else {
    ols_out <- ols(y, z)
    xi <- ols_out$resid
    xi_mat <- matrix(xi, ncol = I)
    n <- nrow(xi_mat)
    Sigma_hat <- (crossprod(xi_mat) / n) %x% diag(n)
    Inv_sigm<-inv_matrix_sqrt(Sigma_hat)
    y_star <- Inv_sigm %*% y
    z_star <- Inv_sigm %*% z
    x_star <- Inv_sigm %*% x

    x1_star <- Inv_sigm %*% x[,1]
    x2_star <- Inv_sigm %*% x[,2]
    x3_star <- Inv_sigm %*% x[,3]

    ols2_out <- ols(y_star, z_star)
    xi_star <- ols2_out$resid

    dayX_1<-ols(x1_star, z_star)
    dayRes<-dayX_1$resid
    weekX_2<-ols(x2_star, z_star)
    weekRes<-weekX_2$resid
    monthX_3<-ols(x3_star, z_star)
    monthRes<-monthX_3$resid

    pi_1<-xi_star*dayRes
    pi_2<-xi_star*weekRes
    pi_3<-xi_star*monthRes

    depvar<-c(rep(1,length(y)))
    Robust_reg<-ols_no_int(depvar,cbind(pi_1,pi_2,pi_3))
    rss<-sum((Robust_reg$fit - depvar) ^ 2)

    ols2_out <- ols(xi_star, cbind(x_star, z_star))
    nu_star <- ols2_out$resid
    # standard Stats
    LM_stat1 <- crossprod(xi_star) - crossprod(nu_star)
    p_value1 <- pchisq(LM_stat1, ncol(x_star), lower.tail = FALSE)
    LM_stat2 <- ((nrow(y_star) - ncol(x_star) - ncol(z_star)) / ncol(x_star)) * LM_stat1 / (nrow(y_star) - LM_stat1)
    p_value2 <- pf(LM_stat2, ncol(x_star), nrow(y_star) - ncol(x_star) - ncol(z_star),
                   lower.tail = FALSE)
    #Heteroscedastic robust stats
    LM_stat3<-(length(y)-rss)
    p_value3 <- pchisq(LM_stat3, ncol(x_star), lower.tail = FALSE)


  }
  LM_out <- matrix(c(LM_stat1, p_value1, LM_stat2, p_value2, LM_stat3, p_value3), nrow = 2, ncol = 3)
  rownames(LM_out) <- c("LM_stat", "p_value")
  colnames(LM_out) <- c("Asymp", "FS_cor", "Asymp_Robust")
  return(LM_out)
}

#' @title Wrapper around APFr
#' @description Wrapper around APFr \code{ apf_plot} to get the desired threshold on the p values such that the FDR is at most \code{ FDR_max}
#' and APF is closest to \code{ APF_min}
#' @param mat_input the input matrix
#' @param FDR_max the maximum false discovery rate desired
#' @param APF_min the minimum average power function desired
#' @param verbose logical for verbose console
#' @return the desired constrained opt result
apf_fdrOpt<-function(mat_input,FDR_max,APF_min,verbose=F){
  if(is.null(FDR_max)){
    index2<-which.min(abs(mat_input[,4]-APF_min))
    result<-mat_input$Gamma[index2]
  }
  else if(is.null(APF_min)){
    index1<-which.min(abs(mat_input[,3]-FDR_max))
    result<-mat_input$Gamma[index1]
  }
  else{
    index1<-which.min(abs(mat_input[,3]-FDR_max))
    index2<-which.min(abs(mat_input[,4]-APF_min))
    if(index1==index2){
      result<-mat_input$Gamma[index1]
    }
    if(index1!=index2){
      index3<-ceiling((index1+index2)/2)
      result<-mat_input$Gamma[index3]
    }
  }
  if(verbose==T){
    Vresult<-list(result,paste("FDR=",mat_input[,3][index1]),paste("APF=",mat_input[,4][index2]))
    return(Vresult)
  }
  if(verbose==F){
    return(result)
  }
}

#'@title simplify output
#'@description This function transforms and simplifies output of the GC testing functions below.
#'@param test_list the test list from \code{ \link{LM_test}}.
#'@param GCpairs the pairs
#'@return the statistics and pvalues
simplify_list <- function(test_list, GCpairs) {
  n_tests <- length(test_list)
  tests <- array(dim = c(2, 2, n_tests))
  selections <- vector("list", length = n_tests)
  for (i in 1:n_tests) {
    tests[, , i] <- test_list[[i]]$tests
    selections[[i]] <- test_list[[i]]$selections
  }
  GC_names <- rep(NA, length(GCpairs))
  for (i in 1:length(GCpairs)) {
    GC_names[i] <- paste(paste(GCpairs[[i]]$GCfrom, collapse = ", "), "->",
                         paste(GCpairs[[i]]$GCto, collapse = ", "))
  }
  dimnames(tests) <- list(stat = c("LM_stat", "p_value"), type = c("Asymp", "FS_cor"),
                          GCtests = GC_names)
  names(selections) <- GC_names
  return(list(tests = tests, selections = selections))
}

#'@title simplify output
#'@description This function transforms and simplifies output of the GC testing functions below.
#'@param test_list the test list from  \code{ \link{LM_test_robust}}.
#'@param GCpairs the pairs
#'@return the statistics and pvalues
simplify_list_RV <- function(test_list, GCpairs) {
  n_tests <- length(test_list)
  tests <- array(dim = c(2, 3, n_tests))
  selections <- vector("list", length = n_tests)
  for (i in 1:n_tests) {
    tests[, , i] <- test_list[[i]]$tests
    selections[[i]] <- test_list[[i]]$selections
  }
  GC_names <- rep(NA, length(GCpairs))
  for (i in 1:length(GCpairs)) {
    GC_names[i] <- paste(paste(GCpairs[[i]]$GCfrom, collapse = ", "), "->",
                         paste(GCpairs[[i]]$GCto, collapse = ", "))
  }
  dimnames(tests) <- list(stat = c("LM_stat", "p_value"), type = c("Asymp", "FS_cor", "Asymp_Robust"),
                          GCtests = GC_names)
  names(selections) <- GC_names
  return(list(tests = tests, selections = selections))
}

#######################################################################################################
####################################### MAIN FUNCTIONS ################################################

#' @title Lag length Selection via BIC empirical upper bound
#' @description Selects the lag length p of the VAR using an empirical upper bound: residuals of
#' the diagonalized VAR are used to build the empirical covariance matrix and an approximation of
#' its determinant that uses the matrix trace is employed to be able to select p using Bayesian Information Criterion
#' @param data a dataframe or matrix of the original set of time series forming the VAR
#' @param p_max maximum lag length to consider, default is 10
#' @return  returns the estimated lag length upper bound
#' @export
#' @examples \dontrun{ lags_upbound_BIC(data,p_max=10)}
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

##############################################  STATIONARY ########################################################

#' @title Test Granger causality in High Dimensional Stationary VARs
#'
#' @param  GCpair     a named list with names GCto and GCfrom containing vectors of the relevant GC variables.
#' @param  data       data matrix or something that can be coerced to a matrix.
#' @param  p          lag length of the VAR
#' @param  bound      lower bound on tuning parameter lambda
#' @param  parallel   TRUE for parallel computing
#' @param  n_cores    nr of cores to use in parallel computing, default is all but one
#' @return            LM test statistics, p values: asymptotic and with finite sample correction and Lasso selections are printed to the console
#' @export
#' @importFrom parallel makeCluster clusterSetRNGStream clusterExport clusterEvalQ detectCores parSapply stopCluster
#' @importFrom stats cor
#' @examples \dontrun{GCpair<-list("GCto"="X1", "GCfrom"="X2")}
#' \dontrun{HDGC_VAR_I0(GCpair,data,p=2,d=2,parallel = T)}
HDGC_VAR_I0 <- function(GCpair, data, p = 1, bound = 0.5 * nrow(data),
                        parallel = FALSE, n_cores = NULL) {
  GCto <- GCpair$GCto #Granger-caused variable
  GCfrom <- GCpair$GCfrom #Granger-causing
  data <- scale(as.matrix(data), scale = FALSE) #scale dataset
  Mat_corr<-as.matrix(cor(data))
  upper<-as.vector(Mat_corr[upper.tri(Mat_corr)])
  amount<-sum(upper >(1-0.001))
  if(any(1-0.001<upper)){
    warning(paste("The used dataset contains",amount, "correlations larger than 0.999, this can cause failure of OLS and poor variable selection.",sep=" "))
  }
  K <- ncol(data) #numb of variables
  X_all <- create_lags(data, p , include.original = FALSE)  #create p lags + d augmentation of all (original K not included) #correspond to lDatafin=xcont
  if((ncol(X_all)+ncol(data))>nrow(data)){
    warning( paste("You are estimating an HD model in which each equation has p*K=",(ncol(X_all)+ncol(data)), " parameters and ", ncol(data), " observations.
                   Depending on how large is p, to avoid failure of OLS you might want to increase the bound=0.5*nrow(data)." ))
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
  lasso_Sy <- active_set_1(i = 1, d=0,p=p, X_index=NULL , y = y_I, z = Z, z_a = NULL, bound = bound)

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

#' @title Test multiple combinations Granger causality in High Dimensional Stationary VARs
#' @description This function is a wrapper around \code{\link{HDGC_VAR_I0}} that allows for multiple combinations to be tested
#' @param data the data matrix or something that can be coerced to a matrix.
#' @param GCpairs it should contain a nested list. The outer list is all the pairs to be considered.
#' The inner list contains the GCto and GCfrom vectors needed for \code{\link{HDGC_VAR_I0}}.
#' @param  p          lag length of the VAR
#' @param  bound      lower bound on tuning parameter lambda
#' @param  parallel   TRUE for parallel computing
#' @param  n_cores    nr of cores to use in parallel computing, default is all but one
#' @return            LM test statistics, p-values: asymptotic and with finite sample correction and Lasso selections are printed to the console
#' @export
#' @importFrom parallel makeCluster clusterSetRNGStream clusterExport clusterEvalQ detectCores parSapply stopCluster parLapply
#' @examples \dontrun{GCpairs<-list(list("GCto"="X1", "GCfrom"="X2"),list("GCto"="X2", "GCfrom"="X3"))}
#' \dontrun{HDGC_VAR_multiple_I0(data, GCpairs, p=1, parallel = T )}
HDGC_VAR_multiple_I0 <- function(data, GCpairs, p = 1, bound = 0.5 * nrow(data),
                              parallel = FALSE, n_cores = NULL) {
  if (parallel) {
    if (is.null(n_cores)) {
      n_cores <- detectCores() - 1
    }
    cl <- makeCluster(n_cores, setup_strategy = "sequential")
    clusterSetRNGStream(cl, sample.int(2^20, size = 1))
    clusterExport(cl = cl, Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv)))
    clusterEvalQ(cl = cl, library(glmnet))

    test_list <- parLapply(cl, GCpairs, HDGC_VAR_I0, data = data, p = p, bound = bound,
                           parallel = FALSE)
    stopCluster(cl)
  } else {
    test_list <- lapply(GCpairs, HDGC_VAR_I0, data = data, p = p, bound = bound,
                        parallel = FALSE)
  }
  out <- simplify_list(test_list, GCpairs)
  return(out)
}

#' @title Test multiple pairs Granger causality in High Dimensional Stationary VARs
#' @description A wrapper around \code{\link{HDGC_VAR_multiple_I0}}. If GCpairs is used, the function is the same as \code{\link{HDGC_VAR_multiple_I0}}.
#' Alternatively, if we want to test all combinations between variables in GCto and GCfrom, these can be given directly. See Example.
#' @param data        the data matrix or something that can be coerced to a matrix.
#' @param  GCpairs it should contain a nested list. The outer list is all the pairs to be considered. See \code{ \link{HDGC_VAR_multiple_I0}}.
#' The inner list contains the GCto and GCfrom vectors needed for \code{\link{HDGC_VAR_I0}}.
#' @param  GCto       all combination variables Granger caused
#' @param  GCfrom     all combination variables Granger causing
#' @param  p          lag length of the VAR
#' @param  bound      lower bound on tuning parameter lambda
#' @param  parallel   TRUE for parallel computing
#' @param  n_cores    nr of cores to use in parallel computing, default is all but one
#' @return            Granger causality matrix and Lasso selections are printed to the console
#' @export
#' @examples \dontrun{GCto = list(c("V1", "V2")); GCfrom = list(c("V3", "V4", "V5"))}
#' \dontrun{HDGC_VAR_multiple_pairs_I0(data,GCto,GCfrom,p=2,parallel = F)}
HDGC_VAR_multiple_pairs_I0 <- function(data, GCpairs = NULL, GCto = NULL, GCfrom = NULL,
                                    p = 1, bound = 0.5 * nrow(data),
                                    parallel = FALSE, n_cores = NULL) {
  varnames <- colnames(data)
  K <- ncol(data)

  direct_to_from <- !(is.null(GCto) | is.null(GCfrom))
  if (direct_to_from) {
    if (!is.null(GCpairs)) {
      warning("Conflicting input given. GCpairs takes precedence")
    } else {
      nGCto <- length(GCto)
      nGCfrom <- length(GCfrom)
      GCpairs <- vector("list", length = nGCto * nGCfrom)
      for (i in 1:nGCto) {
        for (j in 1:nGCfrom) {
          GCpairs[[(i - 1) * nGCfrom + j]] <- list(GCto = GCto[[i]], GCfrom = GCfrom[[j]])
        }
      }
    }
  }
  GC_all_pairs <- HDGC_VAR_multiple_I0(data = data, GCpairs = GCpairs, p = p, bound = bound,
                                    parallel = parallel, n_cores = n_cores)

  GC_matrix <- array(dim = c(K, K, 2, 2))
  dimnames(GC_matrix) <- list(GCto = varnames, GCfrom = varnames,
                              stat = c("LM_stat", "p_value"), type = c("Asymp", "FS_cor"))
  for (i in 1:length(GCpairs)) {
    ind_to <- which(varnames %in% GCpairs[[i]]$GCto)
    ind_from <- which(varnames %in% GCpairs[[i]]$GCfrom)
    GC_matrix[ind_to, ind_from, , ] <- GC_all_pairs$tests[, , i]
  }
  return(list(tests = GC_matrix, selections = GC_all_pairs$selections))
}

#' @title Granger Causality Network in High Dimensional Stationary VARs
#'
#' @description Wrapper around \code{\link{HDGC_VAR_multiple_I0}} which tests Granger causality from each variable to all other variables,
#' one by one. Can therefore be used to construct a network.
#' @param data        the data matrix or something that can be coerced to a matrix.
#' @param  p          lag length of VAR
#' @param  bound      lower bound on tuning parameter lambda
#' @param  parallel   TRUE for parallel computing
#' @param  n_cores    nr of cores to use in parallel computing, default is all but one
#'
#' @return            Granger causality matrix and Lasso selections are printed to the console
#' @export
#' @examples \dontrun{HDGC_VAR_all_I0(data,p=2,parallel=T)}
HDGC_VAR_all_I0 <- function(data, p = 1,  bound = 0.5 * nrow(data),
                         parallel = FALSE, n_cores = NULL) {
  varnames <- colnames(data)
  K <- ncol(data)
  GCpairs <- vector("list", length = K * (K - 1))
  ind <- 0
  for (i in 1:K) {
    for (j in (1:K)[-i]) {
      ind <- ind + 1
      GCpairs[[ind]] <- list(GCto = varnames[i], GCfrom = varnames[j])
    }
  }
  GC_all_pairs <- HDGC_VAR_multiple_I0(data = data, GCpairs = GCpairs, p = p, bound = bound,
                                    parallel = parallel, n_cores = n_cores)
  GC_matrix <- array(dim = c(K, K, 2, 2))
  dimnames(GC_matrix) <- list(GCto = varnames, GCfrom = varnames,
                              stat = c("LM_stat", "p_value"), type = c("Asymp", "FS_cor"))
  for (i in 1:length(GCpairs)) {
    ind_to <- which(varnames %in% GCpairs[[i]]$GCto)
    ind_from <- which(varnames %in% GCpairs[[i]]$GCfrom)
    GC_matrix[ind_to, ind_from, , ] <- GC_all_pairs$tests[, , i]
  }
  return(list(tests = GC_matrix, selections = GC_all_pairs$selections))
}

##############################################  NON STATIONARY ####################################################

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


#' @title Test multiple combinations Granger causality in High Dimensional mixed Integrated and Cointegrated VARs
#' @description This function is a wrapper around \code{\link{HDGC_VAR}} that allows for multiple combinations to be tested
#' @param data the data matrix or something that can be coerced to a matrix.
#' @param GCpairs it should contain a nested list. The outer list is all the pairs to be considered.
#' The inner list contains the GCto and GCfrom vectors needed for \code{\link{HDGC_VAR}}.
#' @param  p          lag length of the VAR
#' @param  d          order of lag augmentation corresponding to suspected max order of integration
#' @param  bound      lower bound on tuning parameter lambda
#' @param  parallel   TRUE for parallel computing
#' @param  n_cores    nr of cores to use in parallel computing, default is all but one
#' @return            LM test statistics, p-values: asymptotic and with finite sample correction and Lasso selections are printed to the console
#' @export
#' @importFrom parallel makeCluster clusterSetRNGStream clusterExport clusterEvalQ detectCores parSapply stopCluster parLapply
#' @examples \dontrun{GCpairs<-list(list("GCto"="X1", "GCfrom"="X2"),list("GCto"="X2", "GCfrom"="X3"))}
#' \dontrun{HDGC_VAR_multiple(data, GCpairs, p=1, d=2, parallel=T)}
HDGC_VAR_multiple <- function(data, GCpairs, p = 1, d = 0, bound = 0.5 * nrow(data),
                              parallel = FALSE, n_cores = NULL) {
  if (parallel) {
    if (is.null(n_cores)) {
      n_cores <- detectCores() - 1
    }
    cl <- makeCluster(n_cores, setup_strategy = "sequential")
    clusterSetRNGStream(cl, sample.int(2^20, size = 1))
    clusterExport(cl = cl, Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv)))
    clusterEvalQ(cl = cl, library(glmnet))

    test_list <- parLapply(cl, GCpairs, HDGC_VAR, data = data, p = p, d = d, bound = bound,
                           parallel = FALSE)
    stopCluster(cl)
  } else {
    test_list <- lapply(GCpairs, HDGC_VAR, data = data, p = p, d = d, bound = bound,
                        parallel = FALSE)
  }
  out <- simplify_list(test_list, GCpairs)
  return(out)
}



#' @title Test multiple pairs Granger causality in High Dimensional mixed Integrated and Cointegrated VARs
#' @description A wrapper around \code{\link{HDGC_VAR_multiple}}. If GCpairs is used,
#' the function is the same as \code{\link{HDGC_VAR_multiple}}. Alternatively, if a we want to test all combinations between
#' variables in GCto and GCfrom, these can be given directly. See Example.
#' @param  data       the data matrix or something that can be coerced to a matrix.
#' @param  GCpairs it should contain a nested list. The outer list is all the pairs to be considered. See \code{ \link{HDGC_VAR_multiple}}.
#' The inner list contains the GCto and GCfrom vectors needed for \code{\link{HDGC_VAR}}.
#' @param  GCto       all combination variables Granger caused
#' @param  GCfrom     all combination variables Granger causing
#' @param  p          lag length of the VAR
#' @param  d          order of lag augmentation corresponding to suspected max order of integration
#' @param  bound      lower bound on tuning parameter lambda
#' @param  parallel   TRUE for parallel computing
#' @param  n_cores    nr of cores to use in parallel computing, default is all but one
#' @return            Granger causality matrix and Lasso selections are printed to the console
#' @export
#' @examples \dontrun{GCto = list(c("V1", "V2")); GCfrom = list(c("V3", "V4", "V5"))}
#' \dontrun{HDGC_VAR_multiple_pairs(data,GCto,GCfrom,p=2,d=2,parallel = T)
#' }
HDGC_VAR_multiple_pairs <- function(data, GCpairs = NULL, GCto = NULL, GCfrom = NULL,
                                    p = 1, d = 0, bound = 0.5 * nrow(data),
                                    parallel = FALSE, n_cores = NULL) {
  varnames <- colnames(data)
  K <- ncol(data)

  direct_to_from <- !(is.null(GCto) | is.null(GCfrom))
  if (direct_to_from) {
    if (!is.null(GCpairs)) {
      warning("Conflicting input given. GCpairs takes precedence")
    } else {
      nGCto <- length(GCto)
      nGCfrom <- length(GCfrom)
      GCpairs <- vector("list", length = nGCto * nGCfrom)
      for (i in 1:nGCto) {
        for (j in 1:nGCfrom) {
          GCpairs[[(i - 1) * nGCfrom + j]] <- list(GCto = GCto[[i]], GCfrom = GCfrom[[j]])
        }
      }
    }
  }
  GC_all_pairs <- HDGC_VAR_multiple(data = data, GCpairs = GCpairs, p = p, d = d, bound = bound,
                                    parallel = parallel, n_cores = n_cores)

  GC_matrix <- array(dim = c(K, K, 2, 2))
  dimnames(GC_matrix) <- list(GCto = varnames, GCfrom = varnames,
                              stat = c("LM_stat", "p_value"), type = c("Asymp", "FS_cor"))
  for (i in 1:length(GCpairs)) {
    ind_to <- which(varnames %in% GCpairs[[i]]$GCto)
    ind_from <- which(varnames %in% GCpairs[[i]]$GCfrom)
    GC_matrix[ind_to, ind_from, , ] <- GC_all_pairs$tests[, , i]
  }
  return(list(tests = GC_matrix, selections = GC_all_pairs$selections))
}



#' @title Granger Causality Network in High Dimensional mixed Integrated and Cointegrated VARs
#'
#' @description Wrapper around \code{\link{HDGC_VAR_multiple}} which tests Granger causality from each variable to all other variables,
#' one by one. Can therefore be used to construct a network.
#' @param data        the data matrix or something that can be coerced to a matrix.
#' @param  p          lag length of VAR
#' @param  d          order of lag augmentation corresponding to suspected max order of integration
#' @param  bound      lower bound on tuning parameter lambda
#' @param  parallel   TRUE for parallel computing
#' @param  n_cores    nr of cores to use in parallel computing, default is all but one
#'
#' @return            Granger causality matrix and Lasso selections are printed to the console
#' @export
#' @examples \dontrun{HDGC_VAR_all(data, p=2, d=2,parallel = T ) }
HDGC_VAR_all <- function(data, p = 1, d = 0, bound = 0.5 * nrow(data),
                         parallel = FALSE, n_cores = NULL) {
  varnames <- colnames(data)
  K <- ncol(data)
  GCpairs <- vector("list", length = K * (K - 1))
  ind <- 0
  for (i in 1:K) {
    for (j in (1:K)[-i]) {
      ind <- ind + 1
      GCpairs[[ind]] <- list(GCto = varnames[i], GCfrom = varnames[j])
    }
  }
  GC_all_pairs <- HDGC_VAR_multiple(data = data, GCpairs = GCpairs, p = p, d = d, bound = bound,
                                    parallel = parallel, n_cores = n_cores)
  GC_matrix <- array(dim = c(K, K, 2, 2))
  dimnames(GC_matrix) <- list(GCto = varnames, GCfrom = varnames,
                              stat = c("LM_stat", "p_value"), type = c("Asymp", "FS_cor"))
  for (i in 1:length(GCpairs)) {
    ind_to <- which(varnames %in% GCpairs[[i]]$GCto)
    ind_from <- which(varnames %in% GCpairs[[i]]$GCfrom)
    GC_matrix[ind_to, ind_from, , ] <- GC_all_pairs$tests[, , i]
  }
  return(list(tests = GC_matrix, selections = GC_all_pairs$selections))
}

################################                                     #######################################
################################ Functions for REALIZED VOLATILITIES #######################################
############################################################################################################
#' @title  Test Granger causality in High Dimensional Stationary Heterogeneous VARs
#'
#' @param  GCpair     A named list with names GCto and GCfrom containing vectors of the relevant GC variables.
#' @param  data       the data matrix or something that can be transformed to a matrix containing realized volatilities.
#' @param  log        default is TRUE, if the realized volatilities are already log transformed then put =FALSE
#' @param  bound      lower bound on tuning parameter lambda
#' @param  parallel   TRUE for parallel computing
#' @param  n_cores    nr of cores to use in parallel computing, default is all but one
#' @return            LM test statistics and p-values: asymptotic, with finite sample correction and asymptotic with heteroscedasticity correction and Lasso selections are printed to the console
#' @export
#' @importFrom parallel makeCluster clusterSetRNGStream clusterExport clusterEvalQ detectCores parSapply stopCluster
#' @importFrom stats cor
#' @examples \dontrun{GCpair<-list("GCto"="X1", "GCfrom"="X2")}
#' \dontrun{HDGC_HVAR(GCpair, data, log=T,parallel = F)}
HDGC_HVAR <- function(GCpair, data, log = TRUE, bound = 0.5 * nrow(data),
                        parallel = FALSE, n_cores = NULL) {
  p = 3 #impose Three lags (HVAR)
  if(log==T){
    data<-log(data)
  }
  GCto <- GCpair$GCto #Granger-caused variable
  GCfrom <- GCpair$GCfrom #Granger-causing
  data <- (as.matrix(data)) #dataset

  Mat_corr<-as.matrix(cor(data))
  upper<-as.vector(Mat_corr[upper.tri(Mat_corr)])
  amount<-sum(upper >(1-0.001))
  if(any(1-0.001<upper)){
    warning(paste("The used dataset contains",amount, "correlations larger than 0.999, this can cause failure of OLS and poor variable selection.",sep=" "))
  }
  K <- ncol(data) #numb of variables
  X_all <- create_lags_RV(data, include.original = FALSE)  #create 3 lags: daily, weekly, monthly aggregates
  if((ncol(X_all)+ncol(data))>nrow(data)){
    warning( paste("You are estimating an HD model in which each equation has p*K=",(ncol(X_all)+ncol(data)), " parameters and ", ncol(data), " observations.
                   Depending on how large is p, to avoid failure of OLS you might want to increase the bound=0.5*nrow(data)." ))
  }
  X <- X_all[, 1:(K * p)] #original K*p lags of regressors , correspond to lDatafin1
  Y <- data[-(1:22), ] #original variables cut the same burn in
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
  Z <- X_all[, -X_index] #all other variables but Granger-causing, corresponds to xcontnogc
  Z_index<-match( GCto_names, colnames(Z) ) #indeces of GCto_names after removing p lags of GCfrom
  Zx <- diag(I) %x% Z   ## takes off the lags of the grangercausing, corresponds to xcontnogc

  # Regressions for X_GC
  if (parallel) {
    if (is.null(n_cores)) {
      n_cores <- detectCores() - 1
    }
    #cl <- makeCluster(n_cores) #this has issues with R studio
    cl <- makeCluster(n_cores, setup_strategy = "sequential")
    clusterSetRNGStream(cl, sample.int(2^20, size = 1))
    clusterExport(cl = cl, Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv)))
    clusterEvalQ(cl = cl, {
      library(glmnet)
      library(zoo)
    })

    lasso_Sx <- parSapply(cl, seq_len(length(X_index)),
                          active_set_1, d=0,p=p, y = X_GC, X_index =NULL, z = Z, z_a =NULL , bound = bound) #regress each time eliminating only the corresponding granger causing lag
    stopCluster(cl)
  }
  else {
    lasso_Sx <- sapply(seq_len(length(X_index)),
                       active_set_1, d=0,p=p, y = X_GC, X_index =NULL, z = Z, z_a =NULL , bound = bound) #regress each time eliminating only the corresponding granger causing lag
  }

  # Regression for y
  lasso_Sy <- active_set_1(i = 1, d=0,p=p, X_index=NULL , y = y_I, z = Z, z_a = NULL, bound = bound)

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

#' @title Test multiple combinations Granger causality in High Dimensional Stationary HVARs
#'
#' @description       This function is a wrapper around \code{\link{HDGC_HVAR}} that allows for multiple combinations to be tested
#' @param  data       the data matrix or something that can be coerced to a matrix containing realized volatilities.
#' @param  GCpairs    it should contain a nested list: the outer list is all the pairs to be considered,
#'                    the inner list contains the GCto and GCfrom vectors needed for \code{ \link{HDGC_HVAR}}. See Example.
#' @param  log        default is TRUE, if the realized volatilities are already log transformed then put to FALSE
#' @param  bound      lower bound on tuning parameter lambda
#' @param  parallel   TRUE for parallel computing
#' @param  n_cores    numberr of cores to use in parallel computing, default is all but one
#' @return            LM test statistics, p-values: asymptotic and with finite sample correction and Lasso selections are printed to the console
#' @export
#' @importFrom parallel makeCluster clusterSetRNGStream clusterExport clusterEvalQ detectCores parSapply stopCluster parLapply
#' @examples \dontrun{GCpairs<-list(list("GCto"="X1", "GCfrom"="X2"),list("GCto"="X2", "GCfrom"="X3"))}
#' \dontrun{HDGC_HVAR_multiple(data,GCpairs,log=T,parallel = T)
#' }
HDGC_HVAR_multiple <- function(data, GCpairs, log = TRUE, bound = 0.5 * nrow(data),
                              parallel = FALSE, n_cores = NULL) {

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

    test_list <- parLapply(cl, GCpairs, HDGC_HVAR, data = data, log=log, bound = bound,
                           parallel = FALSE)
    stopCluster(cl)
  } else {
    test_list <- lapply(GCpairs, HDGC_HVAR, data = data, log=log, bound = bound,
                        parallel = FALSE)
  }
  out <- simplify_list_RV(test_list, GCpairs)
  return(out)
}

#' @title Test multiple pairs Granger causality in High Dimensional Stationary HVARs
#'
#' @description A wrapper around \code{\link{HDGC_HVAR_multiple}}. If GCpairs is used,
#' the function is the same as \code{\link{HDGC_HVAR_multiple}}. Alternatively, if a we want to test all combinations between
#' variables in GCto and GCfrom, these can be given directly. See Example.
#' @param  data       the data matrix or something that can be coerced to a matrix containing realized volatilities
#' @param  GCpairs     it should contain a nested list. The outer list is all the pairs to be considered.
#' The inner list contains the GCto and GCfrom vectors needed for \code{\link{HDGC_HVAR}}. See \code{ \link{HDGC_HVAR_multiple}}.
#' @param  GCto       all combination variables Granger caused
#' @param  GCfrom     all combination variables Granger causing
#' @param  log        default is TRUE, if the realized volatilities are already log transformed then put  to FALSE
#' @param  bound      lower bound on tuning parameter lambda
#' @param  parallel   TRUE for parallel computing
#' @param  n_cores    nr of cores to use in parallel computing, default is all but one
#' @return            Granger causality matrix and Lasso selections are printed to the console
#' @export
#' @examples \dontrun{GCto = list(c("V1", "V2")); GCfrom = list(c("V3", "V4", "V5"))}
#' \dontrun{HDGC_HVAR_multiple_pairs(data,GCto,GCfrom,p=2,d=2,parallel = F)}
HDGC_HVAR_multiple_pairs <- function(data, GCpairs = NULL, GCto = NULL, GCfrom = NULL, log=TRUE,
                                     bound = 0.5 * nrow(data),
                                    parallel = FALSE, n_cores = NULL) {

  varnames <- colnames(data)
  K <- ncol(data)

  direct_to_from <- !(is.null(GCto) | is.null(GCfrom))
  if (direct_to_from) {
    if (!is.null(GCpairs)) {
      warning("Conflicting input given. GCpairs takes precedence")
    } else {
      nGCto <- length(GCto)
      nGCfrom <- length(GCfrom)
      GCpairs <- vector("list", length = nGCto * nGCfrom)
      for (i in 1:nGCto) {
        for (j in 1:nGCfrom) {
          GCpairs[[(i - 1) * nGCfrom + j]] <- list(GCto = GCto[[i]], GCfrom = GCfrom[[j]])
        }
      }
    }
  }
  GC_all_pairs <- HDGC_HVAR_multiple(data = data, GCpairs = GCpairs, log=log, bound = bound,
                                    parallel = parallel, n_cores = n_cores)

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

#' @title Granger Causality Network in High Dimensional Stationary HVARs
#'
#' @description Wrapper around \code{\link{HDGC_HVAR_multiple}} which tests Granger causality from each variable to all other variables,
#' one by one. Can therefore be used to construct a network.
#' @param data        the data matrix or something that can be coerced to a matrix containing realized volatilities
#' @param  log        default is TRUE, if the realized volatilities are already log transformed then put  to FALSE
#' @param  bound      lower bound on tuning parameter lambda
#' @param  parallel   TRUE for parallel computing
#' @param  n_cores    nr of cores to use in parallel computing, default is all but one
#'
#' @return            Granger causality matrix and Lasso selections are printed to the console
#' @export
#' @examples \dontrun{HDGC_HVAR_all(data, log=T, parallel = T) }
HDGC_HVAR_all <- function(data, log = TRUE, bound = 0.5 * nrow(data),
                         parallel = FALSE, n_cores = NULL) {

  varnames <- colnames(data)
  K <- ncol(data)
  GCpairs <- vector("list", length = K * (K - 1))
  ind <- 0
  for (i in 1:K) {
    for (j in (1:K)[-i]) {
      ind <- ind + 1
      GCpairs[[ind]] <- list(GCto = varnames[i], GCfrom = varnames[j])
    }
  }
  GC_all_pairs <- HDGC_HVAR_multiple(data = data, GCpairs = GCpairs, log=log, bound = bound,
                                    parallel = parallel, n_cores = n_cores)
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


###################################################################################################################
############################### REALIZED VOLATILITIES and REALIZED CORRELATIONS ###################################

#' @title Networks of Realized Volatilities conditional on the set of Realized Correlations
#'
#' @param  realized_variances        Dataset of realized volatilities. A matrix or something that can be coerced to a matrix.
#' @param  realized_covariances      Dataset of realized covariances. A matrix or something that can be coerced to a matrix. Note: the columns should exactly
#' be (((ncol(realized_volatilities)^2)-ncol(realized_volatilities))/2)
#' @param  fisher_transf             Logical: if TRUE the correlations are computed and Fisher transformed
#' @param  log                       Default is TRUE, if the realized volatilities are already log transformed then put to FALSE
#' @param  bound                     Lower bound on lambda
#' @param  parallel                  TRUE for parallel computing
#' @param  n_cores                   Nr of cores to use in parallel computing, default is all but one
#' @return   Granger causality matrix and Lasso selections are printed to the console
#' @export
#' @examples \dontrun{ HDGC_HVAR_RV_RCoV_all(real_var, real_cov, fisher_transf=T, log=T,parallel = T) }
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


#' @title Test multiple combinations Granger causality for realized volatilities in High Dimensional Stationary HVARs
#'
#' @description This function is a wrapper around \code{\link{HDGC_HVAR_RVCOV}} that allows for multiple combinations to be tested
#' @param realized_variances Dataset of realized volatilities. A matrix or something that can be coerced to a matrix. Note: the volatilities must not be in logs.
#' @param realized_correlations Dataset of realized correlations. To compute realized correlations from realized variances and realized covariances use \code{\link{Realized_corr}}
#' @param GCpairs it should contain a nested list. The outer list is all the pairs to be considered. See \code{ \link{HDGC_HVAR_RVCOV}}.
#' The inner list contains the GCto and GCfrom vectors needed for \code{\link{HDGC_HVAR_RVCOV}}.
#' @param  log        default is TRUE, if the realized volatilities are already log transformed then put to FALSE
#' @param  bound      lower bound on tuning parameter lambda
#' @param  parallel   TRUE for parallel computing
#' @param  n_cores    nr of cores to use in parallel computing, default is all but one
#' @return            LM test statistics, p-values (asymptotic and with finite sample correction) and Lasso selections are printed to the console
#' @export
#' @importFrom parallel makeCluster clusterSetRNGStream clusterExport clusterEvalQ detectCores parSapply stopCluster
#' @examples \dontrun{ HDGC_HVAR_multiple_RVCOV(real_var, real_corr, GCpairs, log = T,parallel = T)}
HDGC_HVAR_multiple_RVCOV <- function(realized_variances, realized_correlations, GCpairs, log = TRUE, bound = 0.5 * nrow(realized_variances),
                               parallel = FALSE, n_cores = NULL) {

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

    test_list <- parLapply(cl, GCpairs, HDGC_HVAR_RVCOV, realized_variances=realized_variances,
                           realized_correlations=realized_correlations,  bound = bound, parallel = FALSE)
    stopCluster(cl)
  } else {
    test_list <- lapply(GCpairs, HDGC_HVAR_RVCOV, realized_variances=realized_variances,
                        realized_correlations=realized_correlations, bound = bound,
                        parallel = FALSE)
  }
  out <- simplify_list_RV(test_list, GCpairs)
  return(out)
}
############################################################################################################
############################################################################################################

#' @title Plot High-Dimensional Granger causality Networks
#'
#' @param Comb result of function \code{\link{HDGC_VAR_all_I0}}, \code{\link{HDGC_VAR_multiple_pairs_I0}}, \code{\link{HDGC_VAR_all}}, \code{\link{HDGC_VAR_multiple_pairs}}, \code{\link{HDGC_HVAR_all}} or \code{\link{HDGC_HVAR_multiple_pairs}}
#' @param Stat_type either FS_cor (default), Asymp or Asymp_Robust respectively for F small sample correction, standard chi square test, standard chi square test with heteroscedasticity correction
#' @param alpha the desired probability of type one error, default is 0.01.
#' @param multip_corr A list: first element is logical, if TRUE a multiple testing correction using stats::p.adjust() is used. The second
#'                    element of the list define the p.adjust.method=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#'                    "fdr", "none")). If the second element gets the name "APF_FDR" then APFr::apf_fdr which uses empirical Bayes is called and a third and fourth
#'                    elements in the mutip_corr list are required: gamm=c(,,) requires a min, max and step length values to be set for the threshold on the p_values,
#'                    fdr.apf=c(,) requires one or two values: either (NULL,value) or (value,NULL) if one wants to have specified amount of average power (fdr) no matter fdr (average power).
#'                    If both (value,value) are given, the calculated threshold will find the closest combination to both apf and fdr desired. The last element
#'                    of the list is logical: verbose=TRUE if one wants to know how much apf/fdr the testing has.
#' @param ... all parameters for the network plot: see example and igraph documentation.
#' @param cluster A list: first element is logical, if TRUE a cluster plot using igraph::cluster_edge_betweenness() is plotted.
#'                Other elements are respectively: vertex.size, vertex.label.color,vertex.label.cex, vertex.label.dist, edge.curved (see igraph for details).
#' @return a network compatible with igraph
#' @export
#' @importFrom igraph graph_from_adjacency_matrix as.undirected cluster_edge_betweenness V V<-
#' @importFrom APFr apf_fdr apf_plot
#' @importFrom stats p.adjust
#' @importFrom graphics par
#' @examples \dontrun{Plot_GC_all(Comb, "FS_cor",alpha=0.01,multip_corr=list(F), directed=T, layout.circle}
Plot_GC_all<-function(Comb,Stat_type="FS_cor", alpha=0.01, multip_corr=list(F,"bonferroni",gamm = c(1e-04, 0.1, 0.001),fdr.apf=c(0.05,0.6),verb=F),...,
                      cluster=list(F,10,"black",0.51, 1, 0)){
  if(Stat_type=="FS_cor"){
    #build adjaciency matrix
    input<-as.matrix(Comb[["tests"]][,,2,2])}#p_values
  if(Stat_type=="Asymp"){
    #build adjaciency matrix
    input<-as.matrix(Comb[["tests"]][,,2,1])}#p_values
  if(Stat_type=="Asymp_Robust"){
    #build adjaciency matrix
    input<-as.matrix(Comb[["tests"]][,,2,3])}#p_values
  if(multip_corr[[1]]==F){
    input[input < alpha] <- 1 #put =1 values < alpha
    input[is.na(input)] <- 0 #put =0 the diagonal
    input[input != 1] <- 0 #put =0 values > alpha
    network=graph_from_adjacency_matrix(input, mode='directed',diag=F,add.rownames = TRUE )
    V(network)$label = rownames(input)
  }
  if(multip_corr[[1]]==T){
    if(multip_corr[[2]]=="APF_FDR"){
      input_without_NA<-as.vector(input)
      input_without_NA<-input_without_NA[!is.na(input_without_NA)] #take away NAs
      APF_lst<-apf_fdr(data=input_without_NA, type = "pvl", lobs =300 ,
                       seed = 123, gamm = c(multip_corr[[3]]))
      plotAPF_FDR<-apf_plot(APF_lst, tab = TRUE, APF_inf = 0.5, FDR_sup = 0.1)
      alpha_threshold<-apf_fdrOpt(plotAPF_FDR,FDR_max=multip_corr$fdr.apf[1],APF_min=multip_corr$fdr.apf[2],verbose=F)
      if(multip_corr$verb==T){
        apf_fdrOpt(plotAPF_FDR,FDR_max=multip_corr$fdr.apf[1],APF_min=multip_corr$fdr.apf[2],verbose=T)
      }
      input[input < alpha_threshold] <- 1 #put =1 values < alpha
      input[is.na(input)] <- 0 #put =0 the diagonal
      input[input != 1] <- 0 #put =0 values > alpha
      network=graph_from_adjacency_matrix(input, mode='directed',diag=F,add.rownames = TRUE )
      V(network)$label = rownames(input)
    }
    else if (multip_corr[[2]]!="APF_FDR"){
      adj_input<-p.adjust(as.vector(input),method =multip_corr[[2]]) #adjust p-values for multiple testing
      adj_pval_mat<-matrix(adj_input,nrow =nrow(input) ,ncol =ncol(input),byrow = F) #put them back in a matrix
      adj_pval_mat[adj_pval_mat < alpha] <- 1 #put =1 values < alpha
      adj_pval_mat[is.na(adj_pval_mat)] <- 0 #put =0 the diagonal
      adj_pval_mat[adj_pval_mat != 1] <- 0 #put =0 values > alpha
      rownames(adj_pval_mat)<-rownames(input)
      colnames(adj_pval_mat)<-colnames(input)
      network=graph_from_adjacency_matrix(adj_pval_mat, mode='directed',diag=F,add.rownames = TRUE )
      V(network)$label = rownames(adj_pval_mat)
    }
  }
  par(mfrow=c(1,1))
  plot(network,...)
  if(cluster[[1]]==TRUE){
    ## make graph undirected ##
    net.sym <- as.undirected(network, mode= "collapse")
    ceb <- cluster_edge_betweenness(net.sym,directed=T)
    plot(ceb, net.sym, vertex.size=cluster[[2]],vertex.label.color=cluster[[3]],
         vertex.label.cex=cluster[[4]], vertex.label.dist=cluster[[5]], edge.curved=cluster[[6]] )
    #dendPlot(ceb, mode="hclust")
  }
}















