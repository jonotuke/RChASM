#' calculate Dirichlet probabilities with filtering
#'
#' @param nvector observed read counts
#' @param avector alpha vector to change
#' @param cvector change vector for new karyotype
#' @param log return value in log space
#' @param correction error rate for Y mapping for non-Y karyotypes
#'
#' @returns Dirichlet probabilities
ddirichletultinomial <- function(
  nvector,
  avector,
  cvector,
  log = T,
  correction = 1e-6
) {
  Cvector <- base::c(cvector) %>%
    base::replace(. == 0, correction)

  Acvector <- (Cvector * avector) /
    base::sum(Cvector * avector) *
    base::sum(avector)
  N <- base::sum(nvector)

  dens <- extraDistr::ddirmnom(
    base::matrix(nvector, nrow = 1),
    size = N,
    alpha = base::matrix(Acvector, nrow = 1),
    log = log
  )
  return(dens)
}
