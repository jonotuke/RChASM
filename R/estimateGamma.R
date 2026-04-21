#' Estimate contamination from XX and XY based on observed counts
#'
#' @param nvector observed read counts
#' @param avector vector to change
#' @param c1 X multiplier for karyotype 1
#' @param c2 X multiplier for karyotype 2
#' @param correction error rate for Y mapping for non-Y karyotypes

#'
#' @returns estimated contamination
estimateGamma <- function(nvector, avector, c1, c2, correction) {
  n1 <- nvector[1]
  n2 <- nvector[2]
  N <- base::sum(nvector)
  Cvector1 <- base::c(c1, 1) %>%
    base::replace(. == 0, correction)
  Cvector2 <- base::c(c2, 1) %>%
    base::replace(. == 0, correction)
  Avector1 <- (Cvector1 * avector) /
    base::sum(Cvector1 * avector) *
    base::sum(avector)
  Avector2 <- (Cvector2 * avector) /
    base::sum(Cvector2 * avector) *
    base::sum(avector)
  a01 <- base::sum(Avector1)
  E1 <- N * Avector1 / a01
  a02 <- base::sum(Avector2)
  E2 <- N * Avector2 / a02

  W1 <- (E1 + E2)[1]
  W2 <- (E1 + E2)[2]
  w1 <- W1 / (W1 + W2)
  w2 <- W2 / (W1 + W2)

  gamma1 <- (n1 - E2[1]) / (E1[1] - E2[1])
  gamma2 <- (n2 - E2[2]) / (E1[2] - E2[2])

  return(w1 * gamma1 + w2 * gamma2)
}
