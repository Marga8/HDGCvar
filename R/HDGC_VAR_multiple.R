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
