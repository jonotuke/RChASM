#' A function to print the result of the combined analysis to screen
#'
#' @param inChASM the result of a full ChASM analysis (output from runChASM)
#' @param lines the number of lines to print to screen
#'
#' @returns print output
#'
#' @export
printChASM <- function(inChASM, lines = 20) {
  # A function to print the result of the combined analysis to screen
  # inChASM: the result of a full ChASM analysis (output from runChASM)
  # lines: the number of lines to print to screen
  reqNames <- c(
    'sample',
    'protocol',
    'unusual',
    'flags',
    'autosomal_call',
    'sca_call',
    'C_call',
    'autosomal_total',
    'sca_total',
    'automsomal_maxP',
    'sca_maxP'
  )

  if (all(names(inChASM) %in% reqNames)) {
    KAR <- inChASM
  } else {
    if (!("karyotypes" %in% names(inChASM))) {
      stop("input is not the output from a valid RChASM analysis")
    } else {
      KAR <- inChASM$karyotypes
    }
  }

  KAR %>%
    base::print(n = lines)
}
