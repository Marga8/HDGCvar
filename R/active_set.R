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
