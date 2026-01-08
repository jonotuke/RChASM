#' Adjusts alpha vector for contamination karyotypes
#'
#' @param a vector to change
#' @param c1 X multiplier for karyotype 1
#' @param c2 X multiplier for karyotype 2
#' @param g mixing parameter
#' @param correction error rate for Y mapping for non-Y karyotypes
#'
#' @returns adjusted alpha vector
adjustA <- function(a, c1, c2, g, correction = 1e-6) {
  C1 <- base::c(c1, 1) %>%
    base::replace(. == 0, correction)
  C2 <- base::c(c2, 1) %>%
    base::replace(. == 0, correction)
  A1 <- (C1 * a) / base::sum(C1 * a) * base::sum(a)
  A2 <- (C2 * a) / base::sum(C2 * a) * base::sum(a)

  A <- (g * A1 + (1 - g) * A2) / base::sum(g * A1 + (1 - g) * A2)
  return(A)
}
