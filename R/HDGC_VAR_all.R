#' @title Granger Causality Network in High Dimensional mixed Integrated and Cointegrated VARs
#'
#' @description Wrapper around \code{\link{HDGC_VAR_multiple}} which tests Granger causality from each variable to all other variables,
#' one by one. Can therefore be used to construct a network.
#' @param data        the data matrix or something that can be coerced to a matrix.
#' @param  p          lag length of VAR
#' @param  d          order of lag augmentation corresponding to suspected max order of integration
#' @param  bound      lower bound on tuning parameter lambda
#' @param  parallel   TRUE for parallel computing
#' @param  n_cores    nr of cores to use in parallel computing, default is all but one
#'
#' @return            Granger causality matrix and Lasso selections are printed to the console
#' @export
#' @examples \dontrun{HDGC_VAR_all(data, p=2, d=2,parallel = T ) }
HDGC_VAR_all <- function(data, p = 1, d = 0, bound = 0.5 * nrow(data),
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
  GC_all_pairs <- HDGC_VAR_multiple(data = data, GCpairs = GCpairs, p = p, d = d, bound = bound,
                                    parallel = parallel, n_cores = n_cores)
  GC_matrix <- array(dim = c(K, K, 2, 2))
  dimnames(GC_matrix) <- list(GCto = varnames, GCfrom = varnames,
                              stat = c("LM_stat", "p_value"), type = c("Asymp", "FS_cor"))
  for (i in 1:length(GCpairs)) {
    ind_to <- which(varnames %in% GCpairs[[i]]$GCto)
    ind_from <- which(varnames %in% GCpairs[[i]]$GCfrom)
    GC_matrix[ind_to, ind_from, , ] <- GC_all_pairs$tests[, , i]
  }
  return(list(tests = GC_matrix, selections = GC_all_pairs$selections))
}
