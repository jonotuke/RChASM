utils::globalVariables(
  c(
    "automsomal_maxP",
    "autosomal_call",
    "autosomal_total",
    "sca_call",
    "sca_maxP"
  )
)
#' summary_calls
#'
#' @description
#' Combine both autosomal and sex chromosomal analyses into one output
#'
#'
#' @param calls.combined the result of combining both autosomal and sca analyses
#' @param minTotal minimum total for autosomal or sca analyses
#' @param minPosterior minimum value maxP can take (i.e. minimum posterior probability accepted)
#' @param ignoreUnusual filter our unusual observations?
#' @param printProtocol print protocol column?
#'
#' @returns combined calls
#'
#' @export
summary_calls <- function(
  calls.combined,
  minTotal = 6e4,
  minPosterior = 0.95,
  ignoreUnusual = FALSE,
  printProtocol = FALSE
) {
  # Inputs:

  # Define required column names
  reqNames <- base::c(
    'sample',
    'protocol',
    'unusual',
    'flags',
    'autosomal_call',
    'sca_call',
    'autosomal_total',
    'sca_total',
    'automsomal_maxP',
    'sca_maxP'
  )

  # calls.combined must be a data frame or a tibble
  if (
    !checkmate::checkDataFrame(calls.combined) |
      !checkmate::checkTibble(calls.combined)
  ) {
    base::stop('calls.combined must be a dataframe or a tibble.')
  }
  if (!all(reqNames %in% names(calls.combined))) {
    base::stop(paste0(
      'calls.combined is missing columns with names: ',
      paste0(setdiff(reqNames, names(calls.combined)), collapse = ', ')
    ))
  }

  # check that minTotal is a positive integer
  if (!checkmate::checkInt(minTotal, lower = 1)) {
    base::stop("minTotal must be a positive integer.")
  }

  # check that minPosterior is a positive integer
  if (!checkmate::checkDouble(minPosterior, lower = 0, upper = 1)) {
    base::stop("minPosterior must be between zero and one.")
  }

  filteredCalls <- calls.combined %>%
    dplyr::filter(
      (autosomal_call != 'No Aneuploidy') | (!(sca_call %in% c('XX', 'XY'))),
      (autosomal_total > minTotal) | (sca_total > minTotal),
      (automsomal_maxP >= minPosterior) | (sca_maxP >= minPosterior)
    ) %>%
    dplyr::select(
      -tidyr::contains(ifelse(printProtocol, 'placeHolder', 'protocol'))
    )

  if (ignoreUnusual) {
    filteredCalls <- filteredCalls %>%
      dplyr::filter(!unusual)
  }

  return(filteredCalls)
}
