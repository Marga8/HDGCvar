#'@title simplify output
#'@description This function transforms and simplifies output of the GC testing functions below.
#'@keywords internal
#'@param test_list the test list from  \code{ \link{LM_test_robust}}.
#'@param GCpairs the pairs
#'@return the statistics and pvalues
simplify_list_RV <- function(test_list, GCpairs) {
  n_tests <- length(test_list)
  tests <- array(dim = c(2, 3, n_tests))
  selections <- vector("list", length = n_tests)
  for (i in 1:n_tests) {
    tests[, , i] <- test_list[[i]]$tests
    selections[[i]] <- test_list[[i]]$selections
  }
  GC_names <- rep(NA, length(GCpairs))
  for (i in 1:length(GCpairs)) {
    GC_names[i] <- paste(paste(GCpairs[[i]]$GCfrom, collapse = ", "), "->",
                         paste(GCpairs[[i]]$GCto, collapse = ", "))
  }
  dimnames(tests) <- list(stat = c("LM_stat", "p_value"), type = c("Asymp", "FS_cor", "Asymp_Robust"),
                          GCtests = GC_names)
  names(selections) <- GC_names
  return(list(tests = tests, selections = selections))
}
