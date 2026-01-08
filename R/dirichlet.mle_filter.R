#' calculate Dirichlet parameters with filtering
#'
#' @param y input values
#' @param whichK which parameter to return
#'
#' @returns parameters
dirichlet.mle_filter <- function(y, whichK) {
  x <- y[base::apply(y > 0, MARGIN = 1, all), ] %>%
    stats::na.omit()
  return(sirt::dirichlet.mle(x)$alpha[whichK])
}
