#' make c vector for new karyotypes
#'
#' @param x value to replace 1
#' @param y location to make replacement
#' @param n length of vector
#'
#' @returns c vector
makeC <- function(x, y, n) {
  X <- base::rep(1, n)
  if (y == 0) {
    return(X)
  } else {
    X[y] <- x
    return(X)
  }
}
