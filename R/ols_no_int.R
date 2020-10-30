#'@title Ordinary Least Squares without intercept
#'
#'@keywords internal
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
