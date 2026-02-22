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
#' @param inChASM parameter object
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
    stop_message <- sprintf(
      'Samples %s are not in the data (stopping).\n',
      paste0(IDs, collapse = ', ')
    )
    stop(stop_message)
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
      color = ggplot2::guide_legend(
        override.aes = list(fill = 'darkgrey', stroke = 1)
      )
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
