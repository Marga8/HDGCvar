#' Calculate square root of a matrix
#'
#' @keywords internal
#' @param M the matrix
#' @return the matrix square rooted
matrix_sqrt <- function(M) {
  eva_ve <- eigen(M)
  S <- eva_ve$vectors %*% diag(sqrt(eva_ve$values))
  return(S)
}
