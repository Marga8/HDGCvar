#' @title Test multiple combinations Granger causality for realized volatilities in High Dimensional HVARs
#'
#' @description This function is a wrapper around \code{\link{HDGC_HVAR_RVCOV}} that allows for multiple combinations to be tested
#' @param realized_variances Dataset of (stationary) realized volatilities. A matrix or something that can be coerced to a matrix. Note: the volatilities must not be in logs.
#' @param realized_correlations Dataset of (stationary) realized correlations. To compute realized correlations from realized variances and realized covariances use \code{\link{Realized_corr}}
#' @param GCpairs it should contain a nested list. The outer list is all the pairs to be considered. See \code{ \link{HDGC_HVAR_RVCOV}}.
#' The inner list contains the GCto and GCfrom vectors needed for \code{\link{HDGC_HVAR_RVCOV}}.
#' @param  log        default is TRUE, if the realized volatilities are already log transformed then put to FALSE
#' @param  bound      lower bound on tuning parameter lambda
#' @param  parallel   TRUE for parallel computing
#' @param  n_cores    nr of cores to use in parallel computing, default is all but one
#' @return            LM test statistics, p-values (asymptotic and with finite sample correction) and Lasso selections are printed to the console
#' @export
#' @importFrom parallel makeCluster clusterSetRNGStream clusterExport clusterEvalQ detectCores parSapply stopCluster
#' @examples \dontrun{ HDGC_HVAR_multiple_RVCOV(real_var, real_corr, GCpairs, log = TRUE)}
#' @references Hecq, A., Margaritella, L., Smeekes, S., "Granger Causality Testing in High-Dimensional VARs: a Post-Double-Selection Procedure." arXiv preprint arXiv:1902.10991 (2019).
#' @references  Corsi, Fulvio. "A simple approximate long-memory model of realized volatility." Journal of Financial Econometrics 7.2 (2009): 174-196.
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
