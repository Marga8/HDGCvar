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
