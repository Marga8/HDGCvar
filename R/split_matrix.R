#' @title Split Big Matrix
#'
#' @description Splits a big matrix in sub matrices using Kronecker product
#' @param M the matrix
#' @param r the rows of desired cut.
#' @param c the columns of desired cut.
#' @return the splitted matrix
split_matrix <- function(M, r, c){
  nr <- ceiling(nrow(M)/r)
  nc <- ceiling(ncol(M)/c)
  newM <- matrix(NA, nr*r, nc*c)
  newM[1:nrow(M), 1:ncol(M)] <- M

  div_k <- kronecker(matrix(seq_len(nr*nc), nr, byrow = TRUE), matrix(1, r, c))
  matlist <- split(newM, div_k)
  N <- length(matlist)
  mats <- unlist(matlist)
  dim(mats)<-c(r, c, N)
  return(mats)
}
