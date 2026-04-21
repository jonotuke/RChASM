utils::globalVariables(
  c("flags", "P_C")
)
#' combine data
#'
#' @description
#' Combine both autosomal and sex chromosomal analyses into one output
#'
#'
#' @param calls.auto a list or predefined piped string of sample IDs to look for
#' @param calls.sca the autosomal karyotype calls for the samples
#' @param z.scores the sex chromosomal karyotype calls for the samples
#' @param printMissingIDs print the IDs of individuals not in all 3 input files (or just numbers)?
#'
#' @returns combined data
#'
combineData <- function(
  calls.auto,
  calls.sca,
  z.scores,
  printMissingIDs = FALSE
) {
  # Define required column names
  calls.auto.reqNames <- base::c(
    'total',
    'sample',
    'protocol',
    paste0('chr', 1:22)
  )
  calls.sca.reqNames <- base::c('sample', 'P_call', 'px', 'py', 'total', 'maxP')
  z.scores.reqNames <- base::c('sample', 'protocol', 'Nj', 'Nij', 'Zij', 'chr')

  # calls.auto must be a data frame or a tibble
  if (
    !checkmate::checkDataFrame(calls.auto) | !checkmate::checkTibble(calls.auto)
  ) {
    base::stop('calls.auto must be a dataframe or a tibble.')
  }

  if (!all(calls.auto.reqNames %in% names(calls.auto))) {
    base::stop(paste0(
      'calls.auto is missing columns with names: ',
      paste0(setdiff(calls.auto.reqNames, names(calls.auto)), collapse = ', ')
    ))
  }

  # calls.sca must be a data frame or a tibble
  if (
    !checkmate::checkDataFrame(calls.sca) | !checkmate::checkTibble(calls.sca)
  ) {
    base::stop('calls.sca must be a dataframe or a tibble.')
  }
  if (!all(calls.sca.reqNames %in% names(calls.sca))) {
    base::stop(paste0(
      'calls.sca is missing columns with names: ',
      paste0(setdiff(calls.sca.reqNames, names(calls.sca)), collapse = ', ')
    ))
  }

  # z.scores must be a data frame or a tibble
  if (
    !checkmate::checkDataFrame(z.scores) | !checkmate::checkTibble(z.scores)
  ) {
    base::stop('z.scores must be a dataframe or a tibble.')
  }
  if (!all(z.scores.reqNames %in% names(z.scores))) {
    base::stop(paste0(
      'z.scores is missing columns with names: ',
      paste0(setdiff(z.scores.reqNames, names(z.scores)), collapse = ', ')
    ))
  }

  # printMissingIDs must be logical
  if (!checkmate::checkLogical(printMissingIDs)) {
    base::stop('printMissingIDs must be TRUE or FALSE.')
  }

  # Look at number of missing IDs in inputs
  inauto_notsca <- setdiff(calls.auto$sample, calls.sca$sample)
  if (base::length(inauto_notsca) > 0) {
    if (printMissingIDs) {
      error_message <- inauto_notsca %>%
        base::paste0(collapse = ', ') %>%
        paste0('Following IDs in calls.auto but not calls.sca: ', .)
      base::warning(error_message)
    } else {
      error_message <- inauto_notsca %>%
        base::length() %>%
        paste0('Number of IDs from calls.auto missing in calls.sca: ', .)
      base::warning(error_message)
    }
  }
  inauto_notz <- setdiff(calls.auto$sample, unique(z.scores$sample))
  if (base::length(inauto_notz) > 0) {
    if (printMissingIDs) {
      error_message <- inauto_notz %>%
        base::paste0(collapse = ', ') %>%
        paste0('Following IDs in calls.auto but not z.scores: ', .)
      base::warning(error_message)
    } else {
      error_message <- inauto_notz %>%
        base::length() %>%
        paste0('Number of IDs from calls.auto missing in z.scores: ', .)
      base::warning(error_message)
    }
  }
  insca_notauto <- setdiff(calls.sca$sample, calls.auto$sample)
  if (base::length(insca_notauto) > 0) {
    if (printMissingIDs) {
      error_message <- insca_notauto %>%
        base::paste0(collapse = ', ') %>%
        paste0('Following IDs in calls.sca but not calls.auto: ', .)
      base::warning(error_message)
    } else {
      error_message <- insca_notauto %>%
        base::length() %>%
        paste0('Number of IDs from calls.sca missing in calls.auto: ', .)
      base::warning(error_message)
    }
  }
  insca_notz <- setdiff(calls.sca$sample, unique(z.scores$sample))
  if (base::length(insca_notz) > 0) {
    if (printMissingIDs) {
      error_message <- insca_notz %>%
        base::paste0(collapse = ', ') %>%
        paste0('Following IDs in calls.sca but not z.scores: ', .)
      base::warning(error_message)
    } else {
      error_message <- insca_notz %>%
        base::length() %>%
        paste0('Number of IDs from calls.sca missing in z.scores: ', .)
      base::warning(error_message)
    }
  }
  inz_notauto <- setdiff(unique(z.scores$sample), calls.auto$sample)
  if (base::length(inz_notauto) > 0) {
    if (printMissingIDs) {
      error_message <- inz_notauto %>%
        base::paste0(collapse = ', ') %>%
        paste0('Following IDs in z.scores but not calls.auto: ', .)
      base::warning(error_message)
    } else {
      error_message <- inz_notauto %>%
        base::length() %>%
        paste0('Number of IDs from z.scores missing in calls.auto: ', .)
      base::warning(error_message)
    }
  }
  inz_notsca <- setdiff(unique(z.scores$sample), calls.sca$sample)
  if (base::length(inz_notsca) > 0) {
    if (printMissingIDs) {
      error_message <- inz_notsca %>%
        base::paste0(collapse = ', ') %>%
        paste0('Following IDs in z.scores but not calls.sca: ', .)
      base::warning(error_message)
    } else {
      error_message <- inz_notsca %>%
        base::length() %>%
        paste0('Number of IDs from z.scores missing in calls.sca: ', .)
      base::warning(error_message)
    }
  }

  z.scores.fold <- z.scores %>%
    dplyr::filter(stringr::str_detect(
      chr,
      'chr13|chr18|chr21',
      negate = TRUE
    )) %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      .groups = 'keep',
      Xj = base::sum(Zij^2),
      flags = base::sum(abs(Zij) > 2)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(p = 1 - stats::pchisq(Xj, df = 19), unusual = p < 0.05)

  merged.calls <- calls.auto %>%
    dplyr::select(sample, protocol, total, P_call, maxP) %>%
    dplyr::rename(
      autosomal_total = total,
      autosomal_call = P_call,
      automsomal_maxP = maxP
    ) %>%
    dplyr::left_join(
      calls.sca %>%
        dplyr::select(sample, total, P_call, maxP, P_C) %>%
        dplyr::rename(
          sca_total = total,
          sca_call = P_call,
          sca_maxP = maxP,
          C_call = P_C
        ),
      by = c('sample' = 'sample')
    ) %>%
    dplyr::left_join(z.scores.fold, by = c('sample' = 'sample')) %>%
    dplyr::select(
      sample,
      protocol,
      unusual,
      flags,
      dplyr::contains('_call'),
      dplyr::contains('_total'),
      dplyr::contains('_maxP')
    )

  return(merged.calls)
}
