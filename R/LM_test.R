#'@title LM test
#'@description Perform the LM test of significance of x for y given z.
#'@keywords internal
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
