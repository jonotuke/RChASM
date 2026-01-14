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
makeDirichlet <- function(indat, refType, min_reads = 3e4, max_reads = 1e9) {
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
      ggplot2::ggplot(ggplot2::aes(x = px, y = py, col = cluster)) +
      ggplot2::theme_bw() +
      ggplot2::geom_point() +
      ggplot2::xlab('x-rate') +
      ggplot2::ylab('y-rate') +
      ggplot2::scale_colour_discrete(name = 'Clustered\nGenetic\nSex')
    base::plot(a.ggplot)

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
