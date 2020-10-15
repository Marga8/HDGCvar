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
