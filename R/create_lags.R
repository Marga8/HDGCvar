#' @title Lags creation
#'
#' @description Creates a lag matrix of order p. Can include or exclude original series, and trim the NAs in the
#' start of the sample.
#' @param y vector or matrix to lag
#' @param p lag length
#' @param include.original logical, if TRUE the original series are left inside the output matrix.
#' @param trim logical, if TRUE the initial NAs due to the lag gets trimmed
#' @return the lagged vector or matrix
create_lags <- function(y, p = 1, include.original = TRUE, trim = TRUE) {
  x <- as.matrix(y)
  n <- nrow(x)
  k <- ncol(x)
  lx <- matrix(0, nrow = n, ncol = (p + include.original) * k)
  if (is.null(colnames(x))) {
    c.names <- rep("", k)
  } else {
    c.names <- colnames(x)
  }
  colnames(lx) <- rep(c.names, p + include.original)
  for (i in (1 - include.original):p) {
    cols <- k * (i - 1 + include.original) + 1:k
    lx[(1+i):n, cols] <- x[1:(n-i), ]
    colnames(lx)[cols] <- paste(c.names, " l", i, sep = "")
  }
  return(lx[(1 + trim * p):n, , drop = FALSE])
}
