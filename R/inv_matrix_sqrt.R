#' Inverse square root matrix
#'
#' @param M the matrix
#' @return the inverse matrix square rooted
inv_matrix_sqrt <- function(M) {
  if (NROW(M) > 1) {
    eva_ve <- eigen(M)
    S <- diag(sqrt(1/eva_ve$values)) %*% t(eva_ve$vectors)
  } else {
    S <- 1 / sqrt(M)
  }
  return(S)
}
