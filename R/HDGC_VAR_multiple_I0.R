#' @title Test multiple combinations Granger causality in High Dimensional Stationary VARs
#' @description This function is a wrapper around \code{\link{HDGC_VAR_I0}} that allows for multiple combinations to be tested
#' @param data the data matrix or object that can be coerced to a matrix.
#' @param GCpairs it should contain a nested list. The outer list is all the pairs to be considered.
#' The inner list contains the GCto and GCfrom vectors needed for \code{\link{HDGC_VAR_I0}}.
#' @param  p          lag length of the VAR
#' @param  bound      lower bound on tuning parameter lambda
#' @param  parallel   TRUE for parallel computing
#' @param  n_cores    nr of cores to use in parallel computing, default is all but one
#' @return            LM Chi-square test statistics (asymptotic), LM F-stat with finite sample correction, both with their corresponding p-value.
#' Lasso selections are also printed to the console.
#' @export
#' @importFrom parallel makeCluster clusterSetRNGStream clusterExport clusterEvalQ detectCores parSapply stopCluster parLapply
#' @examples \dontrun{GC=list(list("GCto"="Var 1", "GCfrom"="Var 2"),list("GCto"="Var 2", "GCfrom"="Var 3"))}
#' \dontrun{HDGC_VAR_multiple_I0(sample_dataset_I0, GC, p=1 )}
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
