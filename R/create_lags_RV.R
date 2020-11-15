#' @title Lags creation: Daily, Weekly, Monthly aggregation for Realized Volatilities
#'
#' @keywords internal
#' @description Creates a matrix of order 3 containing Daily, Weekly, Monthly returns for Realized Volatilities.
#' Can include or exclude original series, and trim the NAs in the start of the sample. It should be used for HVARs.
#' @param y vector or matrix to lag
#' @param include.original logical, if TRUE the original series are left inside the output matrix.
#' @param trim logical, if TRUE the initial NAs due to the lag gets trimmed
#' @return the lagged vector or matrix
#' @importFrom zoo rollapply
create_lags_RV <- function(y, include.original = TRUE, trim = TRUE) {
  x <- as.matrix(y)
  n <- nrow(x)
  k <- ncol(x)
  lx <- matrix(0, nrow = n, ncol = 4 * k)
  if (is.null(colnames(x))) {
    c.names <- rep("", k)
  } else {
    c.names <- colnames(x)
  }
  colnames(lx) <- rep(c.names, 4)

  for (i in 0:0) {
    cols <-  1:k
    lx[(1+i):n, cols] <- x[1:(n), ]
    colnames(lx)[cols] <- colnames(x)
  }
  for (i in 1:1) {
    cols <- (k+1):(2*k)
    lx[(1+i):n, cols] <- x[1:(n-i), ]
    colnames(lx)[cols] <- paste(c.names, " Daily", " lag", sep = "")
  }
  for (i in 2:2) {
    cols <- ((2*k)+1):(3*k)
    lx[1:5,cols]<-0
    lx[(6):n, cols] <- c(rollapply(x[1:(n-1), ], width = 5, by = 1, FUN = mean, align = "left"))
    colnames(lx)[cols] <- paste(c.names, " Weekly", " lag", sep = "")
  }
  for (i in 3:3) {
    cols <- ((3*k)+1):(4*k)
    lx[1:22, cols]<-0
    lx[(23):n, cols] <- rollapply(x[1:(n-1), ], width = 22, by = 1, FUN = mean, align = "left")
    colnames(lx)[cols] <- paste(c.names, " Monthly", " lag", sep = "")
  }
  if(include.original==F){
    lx<-lx[,(k+1):(4*k)]
  }
  return(lx[(1 + trim * 22):n, , drop = FALSE])
}
