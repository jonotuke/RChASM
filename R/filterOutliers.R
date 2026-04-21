utils::globalVariables(
  c("Outlier", "Obs.Num")
)
#' filter for outliers
#'
#' @param v vector to filter
#'
#' @returns vector for filtering
#'
filterOutliers <- function(v) {
  # Task: filters for outliers
  # Hidden from user
  # Inputs:
  # - v: vector to filter

  if (!checkmate::checkDouble(v)) {
    base::stop('Input vector must be numeric.')
  }

  V <- EnvStats::rosnerTest(
    v,
    k = base::floor(base::length(v) * 0.05),
    warn = F
  )$all.stats %>%
    dplyr::filter(Outlier) %>%
    dplyr::pull(Obs.Num) %>%
    base::sort()
  v[V] <- NA
  return(v)
}
