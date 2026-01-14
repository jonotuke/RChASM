utils::globalVariables(
  c(
    "Nij",
    "Nj",
    "Xj",
    "Zij",
    "alpha",
    "alpha0",
    "chr",
    "flag",
    "muij",
    "p",
    "sigmaij"
  )
)
#' Calculates Z-scores per chromosome
#'
#' @param indat full input tibble
#' @param refType "auto"="autosomal"; "sca"="sex chromosomal"; "diagnostic"
#' @param min_reads min number of reads per sample
#' @param max_reads max number of reads per sample
#'
#' @returns Z-scores
#'
#' @export
makeZscores <- function(
  indat,
  refType = 'auto',
  min_reads = 3e4,
  max_reads = 1e9
) {
  if (!(refType %in% base::c('auto', 'sca', 'diagnostic'))) {
    base::stop(
      "refType must be either 'auto' (Autosomal), 'sca' (Sex Chromosomal) or diagnostic."
    )
  }

  # Set required names
  if (refType == 'auto') {
    reqNames <- base::c(
      'sample',
      'protocol',
      paste0('chr', 1:22),
      'total',
      paste0('p', 1:22)
    )
  } else if (refType == 'sca') {
    reqNames <- base::c(
      'sample',
      'protocol',
      'auto',
      'X',
      'Y',
      'total',
      paste0('p', base::c('x', 'y', 'z'))
    )
  } else if (refType == 'diagnostic') {
    reqNames <- base::c(
      'sample',
      'protocol',
      paste0('chr', 1:22),
      'total',
      paste0('p', 1:22)
    )
  }

  # indat must be a data frame or a tibble
  if (!checkmate::checkDataFrame(indat) | !checkmate::checkTibble(indat)) {
    base::stop('indat must be a dataframe or a tibble.')
  }
  if (!all(reqNames %in% names(indat))) {
    base::stop(base::paste0(
      'indat is missing columns with names: ',
      base::paste0(setdiff(reqNames, base::names(indat)), collapse = ', ')
    ))
  }

  # check that min_reads is a positive integer
  if (!checkmate::checkInt(min_reads, lower = 1)) {
    base::stop("min_reads must be a positive integer.")
  }

  # check that minSamplesPerProtocol is a positive integer
  if (!checkmate::checkInt(max_reads, lower = min_reads)) {
    base::stop(base::paste0(
      "minSamplesPerProtocol must be a positive integer, greater than ",
      min_reads,
      "."
    ))
  }

  inDirichlet <- makeDirichlet(
    indat = indat,
    refType = refType,
    min_reads = min_reads,
    max_reads = max_reads
  )
  alpha_tib <- inDirichlet %>%
    dplyr::select(-a0) %>%
    tidyr::gather(chr, alpha, 2:base::ncol(.)) %>%
    dplyr::mutate(
      chr = chr %>%
        stringr::str_replace_all('a', 'chr')
    )
  z_out <- indat %>%
    dplyr::select(sample, protocol, dplyr::contains('chr')) %>%
    dplyr::rowwise() %>%
    tidyr::gather(chr, Nij, 3:base::ncol(.)) %>%
    dplyr::left_join(
      alpha_tib,
      by = base::c('protocol' = 'protocol', 'chr' = 'chr')
    ) %>%
    stats::na.omit() %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(Nj = sum(Nij)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(sample) %>%
    dplyr::group_by(sample, protocol) %>%
    dplyr::mutate(alpha0 = base::sum(alpha)) %>%
    dplyr::ungroup() %>%
    dplyr::select(sample, protocol, chr, Nij, Nj, alpha, alpha0) %>%
    dplyr::mutate(
      muij = Nj * alpha / alpha0,
      sigmaij = base::sqrt(
        Nj *
          alpha /
          alpha0 *
          (1 - alpha / alpha0) *
          ((Nj + alpha0) / (1 + alpha0))
      ),
      Zij = (Nij - muij) / sigmaij,
      flag = base::abs(Zij) > 2
    ) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(
      Xj = base::sum(Zij^2),
      p = 1 - stats::pchisq(Xj, df = base::ncol(inDirichlet) - 1),
      flags = base::sum(flag),
      unusual = p < 0.05
    ) %>%
    dplyr::ungroup()

  return(dplyr::ungroup(z_out))
}
