##################################################################################################################
##################################################################################################################

utils::globalVariables(
  c("flags")
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
#' @export
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

##################################################################################################################
##################################################################################################################

utils::globalVariables(
  c(
    "protocol_n",
    "arotocol",
    "px",
    "py",
    "cluster",
    "N",
    "ay",
    "a0",
    "ax",
    "az",
    "corr"
  )
)
#' function to make a Dirichlet prior
#'
#' @param indat full input tibble
#' @param refType "auto"="autosomal"; "sca"="sex chromosomal"; "diagnostic"
#' @param min_reads min number of reads per sample
#' @param max_reads max number of reads per sample
#'
#' @returns Dirichlet prior
#'
#' @export
makeDirichlet <- function(
  indat,
  refType,
  min_reads = 3e4,
  max_reads = 1e9,
  show_plot = TRUE
) {
  # refType must be auto, sca or diagnostic
  # Check that refType is one of the two possible inputs
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
    base::stop(paste0(
      'indat is missing columns with names: ',
      paste0(setdiff(reqNames, names(indat)), collapse = ', ')
    ))
  }

  # check that min_reads is a positive integer
  if (!checkmate::checkInt(min_reads, lower = 1)) {
    base::stop("min_reads must be a positive integer.")
  }

  # check that minSamplesPerProtocol is a positive integer
  if (!checkmate::checkInt(max_reads, lower = min_reads)) {
    base::stop(paste0(
      "minSamplesPerProtocol must be a positive integer, greater than ",
      min_reads,
      "."
    ))
  }

  # check that show_plots is boolean
  if (!checkmate::checkLogical(show_plot)) {
    base::stop("show_plot must be TRUE or FALSE")
  }

  # Extract protocols
  protocols <- indat %>%
    dplyr::pull(protocol) %>%
    base::unique()
  # Find column names of interest
  columnIDs <- colnames(indat) %>%
    stringr::str_detect(pattern = '^p[:digit:]+$|^p[xyz]$') %>%
    base::which()

  # Remove chromosomes 13, 18 and 21 for diagnostics
  if (refType == 'diagnostic') {
    columnIDs <- base::setdiff(
      columnIDs,
      base::which(base::colnames(indat) %in% base::c('p13', 'p18', 'p21'))
    )
  }

  # Get protocol counts
  protocolsN <- indat %>%
    dplyr::filter(total >= min_reads, total <= max_reads) %>%
    dplyr::group_by(protocol) %>%
    dplyr::summarise(n = dplyr::n())

  # Protocols require enough observations to estimte priors....
  if (!all(protocolsN$n > length(columnIDs))) {
    base::cat('Protocols without enough observations:\n')
    protocolsN %>%
      dplyr::filter(n < length(columnIDs)) %>%
      dplyr::mutate(protocol_n = paste0(protocol, ' (n=', n, ')')) %>%
      dplyr::pull(protocol_n) %>%
      base::paste0(collapse = ', ') %>%
      base::paste0('\n') %>%
      base::cat() %>%
      base::stop()
  }

  # auto or sca?
  if (refType == 'auto') {
    out_tibble <- tibble::tibble(protocol = protocols)
    outMat <- base::matrix(
      0,
      nrow = base::length(protocols),
      ncol = base::length(columnIDs),
      dimnames = base::list(NULL, base::names(indat)[columnIDs])
    )

    for (i in 1:base::nrow(out_tibble)) {
      outMat[i, 1:base::ncol(outMat)] <- indat %>%
        dplyr::filter(
          protocol == out_tibble$protocol[i],
          total >= min_reads,
          total <= max_reads
        ) %>%
        dplyr::select(tidyr::all_of(columnIDs)) %>%
        dplyr::mutate_all(filterOutliers) %>%
        dirichlet.mle_filter() %>%
        base::unlist()
    }

    out.dirichlet <- outMat %>%
      base::as.data.frame() %>%
      dplyr::bind_cols(out_tibble, .) %>%
      tibble::as_tibble() %>%
      dplyr::rename_at(
        dplyr::vars(tidyr::contains('p')),
        stringr::str_replace_all,
        pattern = 'p',
        replacement = 'a'
      ) %>%
      dplyr::rename(protocol = arotocol) %>%
      dplyr::mutate(a0 = rowSums(dplyr::select(., tidyr::starts_with('a'))))
    return(out.dirichlet)
  } else if (refType == 'diagnostic') {
    out_tibble <- tibble::tibble(protocol = protocols)
    outMat <- base::matrix(
      0,
      nrow = base::length(protocols),
      ncol = base::length(columnIDs),
      dimnames = base::list(NULL, base::names(indat)[columnIDs])
    )

    for (i in 1:nrow(out_tibble)) {
      outMat[i, 1:ncol(outMat)] <- indat %>%
        dplyr::filter(
          protocol == out_tibble$protocol[i],
          total >= min_reads,
          total <= max_reads
        ) %>%
        dplyr::select(tidyr::all_of(columnIDs)) %>%
        dplyr::mutate_all(filterOutliers) %>%
        dirichlet.mle_filter() %>%
        base::unlist()
    }

    out.dirichlet <- outMat %>%
      as.data.frame() %>%
      dplyr::bind_cols(out_tibble, .) %>%
      tibble::as_tibble() %>%
      dplyr::rename_at(
        dplyr::vars(tidyr::contains('p')),
        stringr::str_replace_all,
        pattern = 'p',
        replacement = 'a'
      ) %>%
      dplyr::rename(protocol = arotocol) %>%
      dplyr::mutate(a0 = rowSums(dplyr::select(., tidyr::starts_with('a'))))
    return(out.dirichlet)
  } else if (refType == 'sca') {
    out_tibble <- base::list(
      protocol = protocols,
      sex = base::c('XX', 'XY')
    ) %>%
      base::expand.grid() %>%
      dplyr::arrange(protocol) %>%
      tibble::as_tibble()

    dat.clustered <- indat %>%
      dplyr::filter(total >= min_reads, total < max_reads) %>%
      dplyr::mutate(
        cluster = base::factor(
          mclust::Mclust(tibble::tibble(px = px, py = py), G = 3)$classification
        )
      ) %>%
      dplyr::group_by(cluster) %>%
      dplyr::mutate(
        my = stats::median(py),
        mx = stats::median(px),
        N = dplyr::n()
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(N != base::min(N)) %>%
      dplyr::mutate(
        cluster = dplyr::case_when(
          (my == base::max(my)) & (mx == base::min(mx)) ~ 'XY',
          (my == base::min(my)) & (mx == base::max(mx)) ~ 'XX',
          T ~ 'Unsure'
        )
      ) %>%
      dplyr::filter(cluster != 'Unsure') %>%
      dplyr::mutate(
        cluster = base::factor(
          mclust::Mclust(tibble::tibble(px = px, py = py), G = 2)$classification
        )
      ) %>%
      dplyr::group_by(cluster) %>%
      dplyr::filter(
        !rstatix::mahalanobis_distance(tibble::tibble(px, py))$is.outlier
      ) %>%
      dplyr::mutate(
        my = stats::median(py),
        mx = stats::median(px),
        N = dplyr::n()
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        cluster = dplyr::case_when(
          (my == base::max(my)) & (mx == base::min(mx)) ~ 'XY',
          (my == base::min(my)) & (mx == base::max(mx)) ~ 'XX',
          T ~ 'Unsure'
        )
      )

    a.ggplot <- dat.clustered %>%
      ggplot2::ggplot(ggplot2::aes(
        x = px,
        y = py,
        shape = base::interaction(cluster, protocol, sep = '/'),
        fill = base::interaction(cluster, protocol, sep = '/')
      )) +
      ggplot2::theme_bw() +
      ggplot2::geom_point(col = 'black') +
      ggplot2::xlab('x-rate') +
      ggplot2::ylab('y-rate') +
      ggplot2::scale_fill_discrete(
        name = 'Clustered\nGenetic\nSex\nby Protocol'
      ) +
      ggplot2::scale_shape_manual(
        name = 'Clustered\nGenetic\nSex\nby Protocol',
        values = rep(21:25, 1e3)
      )

    if (show_plot) {
      base::plot(a.ggplot)
    }

    outMat <- base::matrix(
      0,
      nrow = base::length(protocols) * 2,
      ncol = base::length(columnIDs),
      dimnames = base::list(NULL, base::names(indat)[columnIDs])
    )
    for (i in 1:nrow(out_tibble)) {
      outMat[i, 1:ncol(outMat)] <- dat.clustered %>%
        dplyr::filter(
          protocol == out_tibble$protocol[i],
          cluster == out_tibble$sex[i],
          total >= min_reads,
          total <= max_reads
        ) %>%
        dplyr::select(tidyr::all_of(columnIDs)) %>%
        dirichlet.mle_filter() %>%
        base::unlist()
    }
    out.dirichlet <- outMat %>%
      base::as.data.frame() %>%
      dplyr::bind_cols(out_tibble, .) %>%
      tibble::as_tibble() %>%
      dplyr::rename_at(
        dplyr::vars(tidyr::contains('p')),
        stringr::str_replace_all,
        pattern = 'p',
        replacement = 'a'
      ) %>%
      dplyr::rename(protocol = arotocol) %>%
      dplyr::mutate(
        a0 = rowSums(dplyr::select(., tidyr::starts_with('a'))),
        corr = ay / a0
      ) %>%
      dplyr::group_by(protocol) %>%
      dplyr::summarise(
        ax = sum(ax * base::c(0, 1)),
        ay = sum(ay * base::c(0, 1)),
        az = sum(az * base::c(0, 1)),
        a0 = sum(a0 * base::c(0, 1)),
        correction = sum(corr * base::c(1, 0))
      ) %>%
      dplyr::ungroup() %>%
      return(out.dirichlet)
  } else {
    stop('refType must be "auto" or "sca".')
  }
}

##################################################################################################################
##################################################################################################################

utils::globalVariables(
  c(
    "E",
    "auto_total",
    "chrX",
    "chrY",
    "fold",
    "fold_label",
    "isTarget",
    "nij",
    "sca_total",
    "targetID",
    "unusual",
    "x_label"
  )
)
#' Plots diagnostic plots
#'
#' @param IDs a list or predefined piped string of sample IDs to look for
#' @param calls.auto the autosomal karyotype calls for the samples
#' @param calls.sca the sex chromosomal karyotype calls for the samples
#' @param z.scores the autosomal z-scores for the samples
#' @param inDirichlet.auto the dirichlet parameters for the autosomes
#' @param inDirichlet.sca the dirichlet parameters for the sex chromosomes
#' @param min_reads min number of reads per sample
#' @param max_reads max number of reads per sample
#' @param addLabels add "A", "B" and "C" to the plot panels?
#'
#' @returns plots
#'
#' @export
plot_diagnostic <- function(
  IDs,
  inChASM,
  addLabels = FALSE
) {
  calls.auto <- inChASM$karyotypes.auto
  calls.sca <- inChASM$karyotypes.sca
  z.scores <- inChASM$z.scores
  inDirichlet.auto <- inChASM$dirichlet.auto
  inDirichlet.sca <- inChASM$dirichlet.sca
  min_reads = inChASM$min_reads
  max_reads = inChASM$max_reads

  # Define required column names
  calls.auto.reqNames <- base::c(
    'total',
    'sample',
    'protocol',
    paste0('chr', 1:22)
  )
  calls.sca.reqNames <- base::c('sample', 'P_call', 'px', 'py', 'total', 'maxP')
  z.scores.reqNames <- base::c('sample', 'protocol', 'Nj', 'Nij', 'Zij', 'chr')
  inDirichlet.auto.reqNames <- base::c('protocol', paste0('a', 0:22))
  inDirichlet.sca.reqNames <- base::c(
    'protocol',
    'ax',
    'ay',
    'az',
    'a0',
    'correction'
  )

  # IDs must be either a vector or strings or a string
  if (!checkmate::checkCharacter(IDs)) {
    base::stop('IDs must be a string or a vector of strings.')
  }

  Nobs <- z.scores %>%
    dplyr::filter(stringr::str_detect(sample, paste0(IDs, collapse = '|'))) %>%
    nrow() %>%
    {
      . / 22
    }
  if (Nobs == 0) {
    sprintf(
      'Samples %s are not in the data (stopping).\n',
      paste0(IDs, collapse = ', ')
    ) %>%
      cat()
    stop()
  }

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

  # inDirichlet.auto must be a data frame or a tibble
  if (
    !checkmate::checkDataFrame(inDirichlet.auto) |
      !checkmate::checkTibble(inDirichlet.auto)
  ) {
    base::stop('inDirichlet.auto must be a dataframe or a tibble.')
  }
  if (!all(inDirichlet.auto.reqNames %in% names(inDirichlet.auto))) {
    base::stop(paste0(
      'inDirichlet.auto is missing columns with names: ',
      paste0(
        setdiff(inDirichlet.auto.reqNames, names(inDirichlet.auto)),
        collapse = ', '
      )
    ))
  }

  # inDirichlet.sca must be a data frame or a tibble
  if (
    !checkmate::checkDataFrame(inDirichlet.sca) |
      !checkmate::checkTibble(inDirichlet.sca)
  ) {
    base::stop('inDirichlet.sca must be a dataframe or a tibble.')
  }
  if (!all(inDirichlet.sca.reqNames %in% names(inDirichlet.sca))) {
    base::stop(paste0(
      'inDirichlet.sca is missing columns with names: ',
      paste0(
        setdiff(inDirichlet.sca.reqNames, names(inDirichlet.sca)),
        collapse = ', '
      )
    ))
  }

  # check that min_reads is a positive integer
  if (!checkmate::checkInt(min_reads, lower = 1)) {
    base::stop("min_reads must be a positive integer.")
  }

  # check that minSamplesPerProtocol is a positive integer
  if (!checkmate::checkInt(max_reads, lower = min_reads)) {
    base::stop(base::paste0(
      "max_reads must be a positive integer, greater than min_reads",
      min_reads,
      "."
    ))
  }

  gg1 <- calls.sca %>%
    dplyr::filter(
      total >= min_reads,
      total <= max_reads,
      P_call %in% base::c('XX', 'XY'),
      maxP > 0.99
    ) %>%
    dplyr::select(sample, P_call, px, py) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = px,
      y = py,
      col = P_call,
      label = sample
    )) +
    ggplot2::theme_bw() +
    ggplot2::geom_point(pch = 1) +
    ggplot2::geom_point(
      data = calls.sca %>%
        dplyr::select(sample, P_call, px, py) %>%
        dplyr::filter(stringr::str_detect(sample, paste0(IDs, collapse = '|'))),
      pch = 20,
      colour = 'black'
    ) +
    ggrepel::geom_label_repel(
      data = calls.sca %>%
        dplyr::select(sample, P_call, px, py) %>%
        dplyr::filter(stringr::str_detect(sample, paste0(IDs, collapse = '|'))),
      colour = 'black',
      box.padding = 1
    ) +
    ggsci::scale_colour_startrek(drop = FALSE, guide = 'none') +
    ggplot2::xlab('X-Rate') +
    ggplot2::ylab('Y-Rate')

  ## Second ggplot panel
  protocols <- z.scores %>%
    dplyr::filter(stringr::str_detect(
      sample,
      base::paste0(IDs, collapse = '|')
    )) %>%
    dplyr::pull(protocol) %>%
    base::unique()
  IDs_to_keep <- z.scores %>%
    dplyr::filter(protocol %in% protocols) %>%
    dplyr::pull(sample) %>%
    base::unique() %>%
    stringr::str_subset(
      pattern = base::paste0(IDs, collapse = '|'),
      negate = TRUE
    ) %>%
    base::sample(5)
  target_levels <- stringr::str_subset(
    z.scores$sample,
    base::paste0(IDs, collapse = '|')
  ) %>%
    base::unique() %>%
    base::unique() %>%
    base::c('Reference')
  gg2 <- z.scores %>%
    dplyr::filter(
      (Nj >= min_reads & Nj <= max_reads & flags <= 2) |
        stringr::str_detect(sample, base::paste0(IDs, collapse = '|'))
    ) %>%
    dplyr::filter(protocol %in% protocols) %>%
    dplyr::mutate(
      isTarget = stringr::str_detect(sample, paste0(IDs, collapse = '|')),
      targetID = base::ifelse(isTarget, sample, 'Reference') %>%
        base::factor(levels = target_levels)
    ) %>%
    dplyr::group_by(targetID) %>%
    dplyr::filter((isTarget) | (sample %in% IDs_to_keep)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      chr = base::factor(chr, levels = base::paste0('chr', 1:22))
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = chr, y = Zij)) +
    ggplot2::theme_bw() +
    ggplot2::geom_hline(yintercept = base::c(-2, 2), linetype = 'dashed') +
    ggplot2::geom_hline(yintercept = 0, col = 'black') +
    ggstar::geom_star(
      ggplot2::aes(
        starshape = targetID,
        fill = targetID,
        col = isTarget,
        alpha = isTarget
      ),
      size = 3
    ) +
    ggplot2::scale_colour_manual(
      values = base::c('black', 'red'),
      guide = 'none'
    ) +
    ggplot2::scale_alpha_manual(values = base::c(0.2, 1), guide = 'none') +
    ggplot2::scale_fill_discrete(name = NULL) +
    ggstar::scale_starshape_discrete(name = NULL) +
    ggplot2::ylab('Z-score') +
    ggplot2::xlab(NULL) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
      legend.position = 'bottom'
    )

  exp.sca <- inDirichlet.sca %>%
    dplyr::filter(protocol %in% protocols) %>%
    dplyr::mutate(chrX = ax / a0, chrY = ay / a0) %>%
    dplyr::select(protocol, chrX, chrY) %>%
    tidyr::gather(chr, E, chrX, chrY)

  exp.auto <- inDirichlet.auto %>%
    dplyr::filter(protocol %in% protocols) %>%
    dplyr::mutate(
      chr1 = a1 / a0,
      chr2 = a2 / a0,
      chr3 = a3 / a0,
      chr4 = a4 / a0,
      chr5 = a5 / a0,
      chr6 = a6 / a0,
      chr7 = a7 / a0,
      chr8 = a8 / a0,
      chr9 = a9 / a0,
      chr10 = a10 / a0,
      chr11 = a11 / a0,
      chr12 = a12 / a0,
      chr13 = a13 / a0,
      chr14 = a14 / a0,
      chr15 = a15 / a0,
      chr16 = a16 / a0,
      chr17 = a17 / a0,
      chr18 = a18 / a0,
      chr19 = a19 / a0,
      chr20 = a20 / a0,
      chr21 = a21 / a0,
      chr22 = a22 / a0
    ) %>%
    dplyr::select(protocol, chr1:chr22) %>%
    tidyr::gather(chr, E, chr1:chr22)

  exp.values <- dplyr::bind_rows(exp.auto, exp.sca) %>%
    dplyr::mutate(
      chr = chr %>%
        base::factor(levels = base::c(paste0('chr', base::c(1:22, 'X', 'Y'))))
    )

  a1 <- calls.auto %>%
    dplyr::rename(auto_total = total) %>%
    dplyr::select(sample, protocol, auto_total, chr1:chr22) %>%
    dplyr::filter(protocol %in% protocols)
  a2 <- calls.sca %>%
    dplyr::filter(
      ((total >= min_reads) &
        (total <= max_reads) &
        (P_call %in% base::c('XY') & (maxP >= 0.95))) |
        (stringr::str_detect(sample, base::paste0(IDs, collapse = '|')))
    ) %>%
    dplyr::filter(protocol %in% protocols) %>%
    dplyr::rename(sca_total = total) %>%
    dplyr::select(sample, protocol, P_call, sca_total, X, Y)
  gg3.data <- dplyr::left_join(
    a2,
    a1,
    by = base::c('sample' = 'sample', 'protocol' = 'protocol')
  ) %>%
    dplyr::filter(protocol %in% protocols) %>%
    dplyr::left_join(
      z.scores %>%
        dplyr::select(sample, protocol, unusual) %>%
        base::unique(),
      by = base::c('sample' = 'sample', 'protocol' = 'protocol')
    ) %>%
    dplyr::select(
      sample,
      protocol,
      P_call,
      dplyr::contains('_total'),
      unusual,
      X,
      Y,
      dplyr::contains('chr')
    ) %>%
    dplyr::filter(
      (!unusual) |
        (stringr::str_detect(sample, base::paste0(IDs, collapse = '|')))
    ) %>%
    dplyr::mutate(
      isTarget = stringr::str_detect(sample, base::paste0(IDs, collapse = '|')),
      targetID = base::ifelse(isTarget, sample, 'Reference')
    ) %>%
    tidyr::gather('chr', 'Nij', X:chr22) %>%
    dplyr::mutate(
      chr = base::ifelse(
        stringr::str_detect(chr, 'chr'),
        chr,
        base::paste0('chr', chr)
      ) %>%
        base::factor(
          levels = base::c(base::paste0('chr', base::c(1:22, 'X', 'Y')))
        )
    ) %>%
    dplyr::arrange(sample, chr) %>%
    dplyr::filter(targetID != 'Reference') %>%
    dplyr::left_join(
      exp.values,
      by = base::c('protocol' = 'protocol', 'chr' = 'chr')
    ) %>%
    dplyr::mutate(
      nij = base::ifelse(
        stringr::str_detect(chr, 'X|Y'),
        Nij / sca_total,
        Nij / auto_total
      ),
      fold = nij / E * 100,
      fold_label = dplyr::case_when(
        fold >= 300 ~ 300,
        fold <= 50 ~ 50,
        T ~ base::round(fold)
      )
    ) %>%
    dplyr::mutate(
      x_label = base::paste0(sample, '\n', protocol, '\nn=', sca_total)
    ) %>%
    dplyr::left_join(
      z.scores %>%
        dplyr::select(sample, protocol, chr, Zij),
      by = c('sample' = 'sample', 'protocol' = 'protocol', 'chr' = 'chr')
    ) %>%
    dplyr::mutate(
      fill_colour = ifelse(
        abs(Zij) < 2 | base::is.na(Zij),
        'black',
        'lightgrey'
      ),
      text_colour = ifelse(abs(Zij) < 2 | base::is.na(Zij), 'white', 'black')
    ) %>%
    dplyr::mutate(
      chr = chr %>%
        base::factor(
          levels = base::c(base::paste0('chr', base::c(1:22, 'X', 'Y')))
        )
    )

  gg3 <- gg3.data %>%
    ggplot2::ggplot(ggplot2::aes(
      x = chr,
      y = fold_label,
      col = x_label,
      shape = x_label,
      label = base::round(fold)
    )) +
    ggplot2::theme_bw() +
    ggplot2::geom_hline(yintercept = 100, linetype = 'dashed') +
    ggplot2::geom_point(fill = gg3.data$fill_colour, pch = 21, size = 8) +
    ggplot2::geom_text(col = gg3.data$text_colour, size = 3) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    ggplot2::scale_y_log10(
      breaks = seq(50, 300, by = 50),
      limits = base::c(50, 300)
    ) +
    ggplot2::facet_wrap(ggplot2::vars(sample)) +
    ggplot2::scale_shape_manual(values = base::rep(21:25, 1e4), name = NULL) +
    ggplot2::scale_colour_discrete(name = NULL) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab('Fold-increase\n(Compared to XY)') +
    ggplot2::theme(legend.position = 'bottom') +
    ggplot2::guides(
      color = guide_legend(override.aes = list(fill = 'darkgrey', stroke = 1))
    )

  if (addLabels) {
    gg.final <- cowplot::plot_grid(
      gg1,
      gg2,
      nrow = 1,
      labels = base::c('A', 'B'),
      rel_widths = base::c(0.4, 0.7),
      align = 'hv'
    ) %>%
      cowplot::plot_grid(gg3, nrow = 2, labels = base::c(NA, 'C'))
  } else {
    gg.final <- cowplot::plot_grid(
      gg1,
      gg2,
      nrow = 1,
      rel_widths = base::c(0.4, 0.7),
      align = 'hv'
    ) %>%
      cowplot::plot_grid(gg3, nrow = 2)
  }
  return(gg.final)
}

runChASM <- function(
  rawReadCountsIn,
  minSamplesPerProtocol = 30,
  min_reads = 6e4,
  max_reads = 1e9,
  p_contamination = 0.01,
  show_plot = TRUE,
  printMissingIDs = FALSE
) {
  # rawReadCountsIn: the reads counts for each chromosome
  # minSamplesPerProtocol: minimum number of reads per protocol for parameter estimation
  # min_reads: the minimum number of reads for Dirichlet parameter estimation
  # min_reads: the maximum number of reads for Dirichlet parameter estimation
  # p_contamination: the probability of a sample yielding significant contamination
  # show_plot: show the clustering plot for sex chromosomal aneuploidy Dirichlet estimation?
  # printMissingIDs: when combining karyotype calls, return names that are missing (or just the number of missing IDs)?

  if (
    !checkmate::checkDataFrame(rawReadCountsIn) |
      !checkmate::checkTibble(rawReadCountsIn)
  ) {
    base::stop("rawReadCountsIn must be a dataframe or a tibble.")
  }
  if (!("protocol" %in% names(rawReadCountsIn))) {
    base::cat(
      "Warning: no protocol set, assuming all data to come from the same sequencing type (called 'default').\n"
    )
    rawReadCountsIn <- rawReadCountsIn %>% dplyr::mutate(protocol = "default")
  }
  reqNames <- base::c(
    "sample",
    "protocol",
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "X",
    "Y"
  )
  currNames <- reqNames %in% base::names(rawReadCountsIn)
  if (!all(currNames)) {
    base::stop(paste0(
      "rawReadCountsIn is missing columns with names: ",
      paste0(setdiff(reqNames, names(rawReadCountsIn)), collapse = ", ")
    ))
  } else {
    rawReadCountsIn <- rawReadCountsIn %>%
      dplyr::select(tidyr::all_of(reqNames))
  }
  if (!checkmate::checkInt(minSamplesPerProtocol, lower = 1)) {
    base::stop("minSamplesPerProtocol must be a positive integer.")
  }
  if (!checkmate::checkInt(min_reads, lower = 1)) {
    base::stop("min_reads must be a positive integer.")
  }
  if (!checkmate::checkInt(max_reads, lower = min_reads)) {
    base::stop("max_reads must be a positive integer greater than min_reads.")
  }
  if (!checkmate::checkDouble(p_contamination, lower = 0, upper = 1)) {
    base::stop("p_contamination must be a number between zero and one.")
  }
  if (!checkmate::checkLogical(show_plot)) {
    base::stop("show_plot must be TRUE or FALSE")
  }
  if (!checkmate::checkLogical(printMissingIDs)) {
    base::stop("printMissingIDs must be TRUE or FALSE")
  }

  readcounts.auto <- RChASM::processReadCounts(
    rawReadCountsIn = rawReadCountsIn,
    refType = 'auto',
    minSamplesPerProtocol = minSamplesPerProtocol
  )
  readcounts.sca <- RChASM::processReadCounts(
    rawReadCountsIn = rawReadCountsIn,
    refType = 'sca',
    minSamplesPerProtocol = minSamplesPerProtocol
  )

  dirichlet.auto <- RChASM::makeDirichlet(
    indat = readcounts.auto,
    refType = 'auto',
    min_reads = min_reads,
    max_reads = max_reads
  )
  dirichlet.sca <- RChASM::makeDirichlet(
    indat = readcounts.sca,
    refType = 'sca',
    min_reads = min_reads,
    max_reads = max_reads,
    show_plot = show_plot
  )

  z.scores.auto <- RChASM::makeZscores(
    indat = readcounts.auto,
    refType = 'auto',
    min_reads = min_reads,
    max_reads = max_reads
  )

  karyotypes.auto <- RChASM::callKaryotypes(
    indat = readcounts.auto,
    inDirichlet = dirichlet.auto,
    p_contamination = p_contamination
  )
  karyotypes.sca <- RChASM::callKaryotypes(
    indat = readcounts.sca,
    inDirichlet = dirichlet.sca,
    p_contamination = p_contamination
  )
  karyotypes <- RChASM::combineData(
    calls.auto = karyotypes.auto,
    calls.sca = karyotypes.sca,
    z.scores = z.scores.auto,
    printMissingIDs = printMissingIDs
  )

  base::
  return(base::list(
    karyotypes = karyotypes,
    karyotypes.auto = karyotypes.auto,
    karyotypes.sca = karyotypes.sca,
    dirichlet.auto = dirichlet.auto,
    dirichlet.sca = dirichlet.sca,
    z.scores = z.scores.auto,
    minSamplesPerProtocol = 30,
    min_reads = 6e4,
    max_reads = 1e9,
    p_contamination = 0.01
  ))
}

##################################################################################################################
##################################################################################################################

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
  inChASM,
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

  calls.combined <- inChASM$karyotypes

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

saveChASM <- function(inChASM, file, sort_by_samplename = FALSE) {
  # A function for saving the results of a ChASM analysis as a tsv
  # inChASM: the result of a full ChASM analysis (output from runChASM)
  # file: the full path and file name to place output
  # sort_by_samplename: reorder alphabetically by sample names?
  if (!checkmate::checkLogical(sort_by_samplename)) {
    base::stop("sort_by_samplename must be TRUE or FALSE")
  }
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

  if (sort_by_samplename) {
    KAR %>%
      dplyr::arrange(sample) %>%
      readr::write_tsv(file = file)
  } else {
    KAR %>%
      readr::write_tsv(file = file)
  }
}
