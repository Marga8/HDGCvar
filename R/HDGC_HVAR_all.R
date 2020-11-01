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
#' @examples \dontrun{HDGC_HVAR_all(data=sample_RV, log=TRUE, parallel = TRUE) }
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
