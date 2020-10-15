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
