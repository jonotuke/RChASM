#' bounds number to be in \[0,1\] (contamination estimation)
#'
#' @param x number
#'
#' @returns bounded number
bound <- function(x) {
  if (!checkmate::checkDouble(x)) {
    base::stop('x must be a number.')
  }
  return(base::min(base::max(x, 0), 1))
}
