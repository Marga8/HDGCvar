#' @title Test multiple pairs Granger causality in High Dimensional HVARs
#'
#' @description A wrapper around \code{\link{HDGC_HVAR_multiple}}. If GCpairs is used,
#' the function is the same as \code{\link{HDGC_HVAR_multiple}}. Alternatively, if a we want to test all combinations between
#' variables in GCto and GCfrom, these can be given directly. See Example.
#' @param  data       the data matrix or something that can be coerced to a matrix containing (stationary) realized volatilities
#' @param  GCpairs     it should contain a nested list. The outer list is all the pairs to be considered.
#' The inner list contains the GCto and GCfrom vectors needed for \code{\link{HDGC_HVAR}}. See \code{ \link{HDGC_HVAR_multiple}}.
#' @param  GCto       all combination variables Granger caused
#' @param  GCfrom     all combination variables Granger causing
#' @param  log        default is TRUE, if the realized volatilities are already log transformed then put  to FALSE
#' @param  bound      lower bound on tuning parameter lambda
#' @param  parallel   TRUE for parallel computing
#' @param  n_cores    nr of cores to use in parallel computing, default is all but one
#' @return            Granger causality matrix and Lasso selections are printed to the console
#' @export
#' @examples \dontrun{GCto = list(c("Var 1", "Var 2")); GCfrom = list(c("Var 3", "Var 4", "Var 5"))}
#' \dontrun{HDGC_HVAR_multiple_pairs(sample_RV,GCto,GCfrom, log=TRUE)}
#' @references Hecq, A., Margaritella, L., Smeekes, S., "Granger Causality Testing in High-Dimensional VARs: a Post-Double-Selection Procedure." arXiv preprint arXiv:1902.10991 (2019).
#' @references  Corsi, Fulvio. "A simple approximate long-memory model of realized volatility." Journal of Financial Econometrics 7.2 (2009): 174-196.
HDGC_HVAR_multiple_pairs <- function(data, GCpairs = NULL, GCto = NULL, GCfrom = NULL, log=TRUE,
                                     bound = 0.5 * nrow(data),
                                     parallel = FALSE, n_cores = NULL) {

  varnames <- colnames(data)
  K <- ncol(data)

  direct_to_from <- !(is.null(GCto) | is.null(GCfrom))
  if (direct_to_from) {
    if (!is.null(GCpairs)) {
      warning("Conflicting input given. GCpairs takes precedence")
    } else {
      nGCto <- length(GCto)
      nGCfrom <- length(GCfrom)
      GCpairs <- vector("list", length = nGCto * nGCfrom)
      for (i in 1:nGCto) {
        for (j in 1:nGCfrom) {
          GCpairs[[(i - 1) * nGCfrom + j]] <- list(GCto = GCto[[i]], GCfrom = GCfrom[[j]])
        }
      }
    }
  }
  GC_all_pairs <- HDGC_HVAR_multiple(data = data, GCpairs = GCpairs, log=log, bound = bound,
                                     parallel = parallel, n_cores = n_cores)

  # GC_matrix <- array(dim = c(K, K, 2, 3))
  # dimnames(GC_matrix) <- list(GCto = varnames, GCfrom = varnames,
  #                             stat = c("LM_stat", "p_value"), type = c("Asymp", "FS_cor", "Asymp_Robust"))
  # for (i in 1:length(GCpairs)) {
  #   ind_to <- which(varnames %in% GCpairs[[i]]$GCto)
  #   ind_from <- which(varnames %in% GCpairs[[i]]$GCfrom)
  #   GC_matrix[ind_to, ind_from, , ] <- GC_all_pairs$tests[, , i]
  # }
  return(list(tests = GC_all_pairs$tests, selections = GC_all_pairs$selections))
}
