############################
# Karyotype functions      -
# Created by: AB Rohrlach  -
# Created on: 10-10-2025   -
############################

#### required packages ----
pacman::p_load(
  tidyverse,
  ggplot2,
  stringr,
  checkmate,
  ggsci,
  ggrepel,
  ggstar,
  rstatix,
  EnvStats,
  extraDistr
)

#### function descriptions ---
processReadCounts <- function(
  rawReadCountsIn,
  refType,
  minTotal = 1e3,
  minSamplesPerProtocol = 30
) {
  # Task: takes read counts and makes useful file
  # Inputs:
  # - rawReadCountsIn: the reads counts via KP
  # - onTarget: for capture, use only on-target (F)?
  # - minSamplesPerProtocol: minimum number of protocol counts to be included
  # - refType: "auto"="autosomal"; "sca"="sex chromosomal"
  # Checks made: yes

  # Need to check inputs
  # Check rawReadCountsIn is a dataframe or tibble
  if (
    !checkmate::checkDataFrame(rawReadCountsIn) |
      !checkmate::checkTibble(rawReadCountsIn)
  ) {
    base::stop('rawReadCountsIn must be a dataframe or a tibble.')
  }

  # Check that rawReadCountsIn contains the minimum required columns.
  # If protocol is missing, add a default
  # Then make sure only required columns are carried through
  if (!('protocol' %in% names(rawReadCountsIn))) {
    base::cat(
      "Warning: no protocol set, assuming all data to come from the same sequencing type (called 'default').\n"
    )
    rawReadCountsIn <- rawReadCountsIn %>%
      dplyr::mutate(protocol = 'default')
  }
  reqNames <- base::c(
    'sample',
    'protocol',
    'chr1',
    'chr2',
    'chr3',
    'chr4',
    'chr5',
    'chr6',
    'chr7',
    'chr8',
    'chr9',
    'chr10',
    'chr11',
    'chr12',
    'chr13',
    'chr14',
    'chr15',
    'chr16',
    'chr17',
    'chr18',
    'chr19',
    'chr20',
    'chr21',
    'chr22',
    'X',
    'Y'
  )
  currNames <- reqNames %in% base::names(rawReadCountsIn)
  if (!all(currNames)) {
    base::stop(paste0(
      'rawReadCountsIn is missing columns with names: ',
      paste0(setdiff(reqNames, names(rawReadCountsIn)), collapse = ', ')
    ))
  } else {
    rawReadCountsIn <- rawReadCountsIn %>%
      dplyr::select(all_of(reqNames))
  }

  # Check that reyType is one of the two possible inputs
  if (!(refType %in% base::c('auto', 'sca'))) {
    base::stop(
      "refType must be either 'auto' (Autosomal) or 'sca' (Sex Chromosomal)."
    )
  }

  # check that minTotal is a positive integer
  if (!checkmate::checkInt(minTotal, lower = 1)) {
    base::stop("minTotal must be a positive integer.")
  }

  # check that minSamplesPerProtocol is a positive integer
  if (!checkmate::checkInt(minSamplesPerProtocol, lower = 1)) {
    base::stop("minSamplesPerProtocol must be a positive integer.")
  }

  # Filter data for minimum values
  rawReadCounts <- rawReadCountsIn %>%
    dplyr::mutate(protocol = tolower(protocol)) %>%
    dplyr::group_by(protocol) %>%
    dplyr::filter(n() >= minSamplesPerProtocol) %>%
    dplyr::ungroup()

  # If refType is SCA then fold into X, Y and Automsomal
  if (tolower(refType) == "sca") {
    outReadCounts <- rawReadCounts %>%
      dplyr::mutate(auto = rowSums(dplyr::select(., starts_with('chr')))) %>%
      dplyr::select(sample, protocol, auto, X, Y) %>%
      dplyr::mutate(total = auto + X + Y) %>%
      dplyr::filter(total >= minTotal) %>%
      dplyr::mutate(px = X / total, py = Y / total, pz = auto / total)
  } else if (tolower(refType) == 'auto') {
    outReadCounts <- rawReadCounts %>%
      dplyr::select(-X, -Y) %>%
      dplyr::mutate(total = rowSums(dplyr::select(., starts_with('chr')))) %>%
      dplyr::filter(total >= minTotal) %>%
      dplyr::select(sample, protocol, total, chr1:chr22) %>%
      dplyr::mutate(
        p1 = chr1 / total,
        p2 = chr2 / total,
        p3 = chr3 / total,
        p4 = chr4 / total,
        p5 = chr5 / total,
        p6 = chr6 / total,
        p7 = chr7 / total,
        p8 = chr8 / total,
        p9 = chr9 / total,
        p10 = chr10 / total,
        p11 = chr11 / total,
        p12 = chr12 / total,
        p13 = chr13 / total,
        p14 = chr14 / total,
        p15 = chr15 / total,
        p16 = chr16 / total,
        p17 = chr17 / total,
        p18 = chr18 / total,
        p19 = chr19 / total,
        p20 = chr20 / total,
        p21 = chr21 / total,
        p22 = chr22 / total
      )
  }
  base::
  return(outReadCounts)
}

makeDirichlet <- function(indat, refType, min_reads = 3e4, max_reads = 1e9) {
  # Task: function to make a Dirichlet prior
  # Inputs:
  # - in_dat: full input tibble
  # - refType: "auto"="autosomal"; "sca"="sex chromosomal"; "diagnostic"
  # - min_reads: min number of reads per sample
  # - max_reads: max number of reads per sample
  # Checks made: yes

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
    dplyr::summarise(n = n())

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
        dplyr::select(all_of(columnIDs)) %>%
        dplyr::mutate_all(filterOutliers) %>%
        dirichlet.mle_filter() %>%
        base::unlist()
    }

    out.dirichlet <- outMat %>%
      base::as.data.frame() %>%
      dplyr::bind_cols(out_tibble, .) %>%
      tibble::as_tibble() %>%
      dplyr::rename_at(
        dplyr::vars(contains('p')),
        str_replace_all,
        pattern = 'p',
        replacement = 'a'
      ) %>%
      dplyr::rename(protocol = arotocol) %>%
      dplyr::mutate(a0 = rowSums(dplyr::select(., starts_with('a'))))
    base::
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
        dplyr::select(all_of(columnIDs)) %>%
        dplyr::mutate_all(filterOutliers) %>%
        dirichlet.mle_filter() %>%
        base::unlist()
    }

    out.dirichlet <- outMat %>%
      as.data.frame() %>%
      dplyr::bind_cols(out_tibble, .) %>%
      as_tibble() %>%
      dplyr::rename_at(
        vars(contains('p')),
        str_replace_all,
        pattern = 'p',
        replacement = 'a'
      ) %>%
      dplyr::rename(protocol = arotocol) %>%
      dplyr::mutate(a0 = rowSums(dplyr::select(., starts_with('a'))))
    base::
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
        dplyr::select(all_of(columnIDs)) %>%
        dirichlet.mle_filter() %>%
        base::unlist()
    }
    out.dirichlet <- outMat %>%
      base::as.data.frame() %>%
      dplyr::bind_cols(out_tibble, .) %>%
      tibble::as_tibble() %>%
      dplyr::rename_at(
        dplyr::vars(contains('p')),
        str_replace_all,
        pattern = 'p',
        replacement = 'a'
      ) %>%
      dplyr::rename(protocol = arotocol) %>%
      dplyr::mutate(
        a0 = rowSums(dplyr::select(., starts_with('a'))),
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

callKaryotypes <- function(indat, inDirichlet, p_contamination = 0.1) {
  # Task: take prior and observed counts and calls karyotypes and Z-scores
  # Inputs:
  # - indat: full input tibble
  # - inDirichlet: associated Dirichlet parameters
  # - p_contamination: the assumed probability of contamination
  # Checks made: yes

  # Is it SCA or Auto?
  refType <- base::ifelse(any(str_detect(names(indat), '^px$')), 'sca', 'auto')

  # indat has to be a data frame or tibble
  if (refType == 'auto') {
    reqNamesIndat <- base::c(
      'sample',
      'protocol',
      base::paste0('chr', 1:22),
      'total',
      paste0('p', 1:22)
    )
    reqNamesDirichlet <- base::c('protocol', paste0('a', 1:22), 'a0')
  } else if (refType == 'sca') {
    reqNamesIndat <- base::c(
      'sample',
      'protocol',
      'auto',
      'X',
      'Y',
      'total',
      paste0('p', base::c('x', 'y', 'z'))
    )
    reqNamesDirichlet <- base::c(
      'protocol',
      'ax',
      'ay',
      'az',
      'a0',
      'correction'
    )
  } else if (refType == 'diag') {
    reqNamesIndat <- base::c(
      'sample',
      'protocol',
      base::paste0('chr', 1:22),
      'total',
      base::paste0('p', 1:22)
    )
    reqNamesDirichlet <- base::c('protocol', paste0('a', 1:22), 'a0')
  }

  # indat must be a data frame or a tibble
  if (!checkmate::checkDataFrame(indat) | !checkmate::checkTibble(indat)) {
    base::stop('indat must be a dataframe or a tibble.')
  }
  if (!all(reqNamesIndat %in% names(indat))) {
    base::stop(paste0(
      'indat is missing columns with names: ',
      paste0(setdiff(reqNamesIndat, names(indat)), collapse = ', ')
    ))
  }

  # inDirichlet must be a data frame or a tibble
  if (
    !checkmate::checkDataFrame(inDirichlet) |
      !checkmate::checkTibble(inDirichlet)
  ) {
    base::stop('inDirichlet must be a dataframe or a tibble.')
  }
  if (!all(reqNamesDirichlet %in% names(inDirichlet))) {
    base::stop(paste0(
      'inDirichlet is missing columns with names: ',
      paste0(setdiff(reqNamesIndat, names(reqNamesDirichlet)), collapse = ', ')
    ))
  }

  # check that p_contamination is a value between 0 and 1
  if (!checkmate::checkDouble(p_contamination, lower = 0, upper = 1)) {
    base::stop("p_contamination must be a number between zero and one.")
  }

  if (refType == 'sca') {
    # SCA calls
    pc <- p_contamination
    p.XXY <- (105 / 205) / 750
    p.X <- (100 / 205) / 2500
    p.XXX <- (100 / 205) / 1000
    p.XYY <- (105 / 205) / 1000
    p.XY <- (1 - p.XXY - p.X - p.XXX - p.XYY) * (105 / 205)
    p.XX <- (1 - p.XXY - p.X - p.XXX - p.XYY) * (100 / 205)

    out_calls <- indat %>%
      dplyr::left_join(
        inDirichlet,
        by = base::c('protocol' = 'protocol'),
        relationship = 'many-to-many'
      ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        Gamma = estimateGamma(
          nvector = base::c(X, Y, auto),
          avector = base::c(ax, ay, az),
          c1 = base::c(1, 1),
          c2 = base::c(2, 0),
          correction = correction
        ) %>%
          bound(),
        P_XY = ddirichletultinomial(
          nvector = base::c(X, Y, auto),
          avector = base::c(ax, ay, az),
          cvector = base::c(1, 1, 1),
          log = T,
          correction = correction
        ) +
          base::log(p.XY),
        P_XX = ddirichletultinomial(
          nvector = base::c(X, Y, auto),
          avector = base::c(ax, ay, az),
          cvector = base::c(2, 0, 1),
          log = T,
          correction = correction
        ) +
          base::log(p.XX),
        P_XXY = ddirichletultinomial(
          nvector = base::c(X, Y, auto),
          avector = base::c(ax, ay, az),
          cvector = base::c(2, 1, 1),
          log = T,
          correction = correction
        ) +
          base::log(p.XXY),
        P_X = ddirichletultinomial(
          nvector = base::c(X, Y, auto),
          avector = base::c(ax, ay, az),
          cvector = base::c(1, 0, 1),
          log = T,
          correction = correction
        ) +
          base::log(p.X),
        P_XXX = ddirichletultinomial(
          nvector = base::c(X, Y, auto),
          avector = base::c(ax, ay, az),
          cvector = base::c(3, 0, 1),
          log = T,
          correction = correction
        ) +
          base::log(p.XXX),
        P_XYY = ddirichletultinomial(
          nvector = base::c(X, Y, auto),
          avector = base::c(ax, ay, az),
          cvector = base::c(1, 2, 1),
          log = T,
          correction = correction
        ) +
          base::log(p.XYY),
        P_cont = ddirichletultinomial(
          nvector = base::c(X, Y, auto),
          avector = adjustA(
            base::c(ax, ay, az),
            base::c(1, 1),
            base::c(2, 0),
            Gamma,
            correction = 1e-6
          ),
          cvector = base::c(1, 1, 1),
          log = T,
          correction = correction
        ) +
          base::log(pc) +
          base::log(p.XY) +
          base::log(p.XX)
      ) %>%
      dplyr::mutate(
        SumP = matrixStats::logSumExp(base::c(
          P_XX,
          P_XY,
          P_XXY,
          P_X,
          P_XXX,
          P_XYY
        )) %>%
          base::exp(),
        SumC = matrixStats::logSumExp(base::c(
          P_XX + base::log(1 - pc),
          P_XY + base::log(1 - pc),
          P_XYY + base::log(1 - pc),
          P_cont
        )) %>%
          base::exp(),
        Posterior_XY = base::exp(P_XY) / SumP,
        Posterior_XX = base::exp(P_XX) / SumP,
        Posterior_XXY = base::exp(P_XXY) / SumP,
        Posterior_X = base::exp(P_X) / SumP,
        Posterior_XXX = base::exp(P_XXX) / SumP,
        Posterior_XYY = base::exp(P_XYY) / SumP
      ) %>%
      dplyr::mutate(
        P_call = dplyr::case_when(
          P_XX ==
            base::max(base::c(P_XX, P_XY, P_XXY, P_X, P_XXX, P_XYY)) ~ 'XX',
          P_XY ==
            base::max(base::c(P_XX, P_XY, P_XXY, P_X, P_XXX, P_XYY)) ~ 'XY',
          P_XXY ==
            base::max(base::c(P_XX, P_XY, P_XXY, P_X, P_XXX, P_XYY)) ~ 'XXY',
          P_X ==
            base::max(base::c(P_XX, P_XY, P_XXY, P_X, P_XXX, P_XYY)) ~ 'X0',
          P_XXX ==
            base::max(base::c(P_XX, P_XY, P_XXY, P_X, P_XXX, P_XYY)) ~ 'XXX',
          P_XYY ==
            base::max(base::c(P_XX, P_XY, P_XXY, P_X, P_XXX, P_XYY)) ~ 'XYY'
        ),
        maxP = base::max(base::c(
          Posterior_XX,
          Posterior_XY,
          Posterior_XXY,
          Posterior_X,
          Posterior_XXX,
          Posterior_XYY
        )),
        Posterior_cont = base::exp(P_cont) / SumC,
        P_C = base::ifelse(
          (Posterior_cont > Posterior_XX) &
            (Posterior_cont > Posterior_XY) &
            (Posterior_cont > Posterior_XXY) &
            (Posterior_cont > Posterior_XYY) &
            (Posterior_cont > Posterior_XXX) &
            (Posterior_cont > Posterior_X) &
            (Gamma != 0) &
            (Gamma != 1),
          'Possible Contamination',
          'No Sig XX+XY'
        )
      )

    final_out_calls <- out_calls %>%
      dplyr::select(
        sample,
        protocol,
        total,
        P_call,
        maxP,
        P_cont,
        dplyr::everything()
      )
    return(final_out_calls)
  } else {
    # Auto calls
    p.13 <- 1 / 7143
    p.18 <- 1 / 3226
    p.21 <- 1 / 705
    p.n <- 1 - base::sum(base::c(p.13, p.18, p.21))

    out_calls <- indat %>%
      dplyr::left_join(inDirichlet, by = base::c('protocol' = 'protocol')) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        P_N = ddirichletultinomial(
          nvector = base::c(
            chr1,
            chr2,
            chr3,
            chr4,
            chr5,
            chr6,
            chr7,
            chr8,
            chr9,
            chr10,
            chr11,
            chr12,
            chr13,
            chr14,
            chr15,
            chr16,
            chr17,
            chr18,
            chr19,
            chr20,
            chr21,
            chr22
          ),
          avector = base::c(
            a1,
            a2,
            a3,
            a4,
            a5,
            a6,
            a7,
            a8,
            a9,
            a10,
            a11,
            a12,
            a13,
            a14,
            a15,
            a16,
            a17,
            a18,
            a19,
            a20,
            a21,
            a22
          ),
          cvector = makeC(1, 0, 22),
          log = T,
          correction = 0
        ) +
          base::log(p.n),
        P_13 = ddirichletultinomial(
          nvector = base::c(
            chr1,
            chr2,
            chr3,
            chr4,
            chr5,
            chr6,
            chr7,
            chr8,
            chr9,
            chr10,
            chr11,
            chr12,
            chr13,
            chr14,
            chr15,
            chr16,
            chr17,
            chr18,
            chr19,
            chr20,
            chr21,
            chr22
          ),
          avector = base::c(
            a1,
            a2,
            a3,
            a4,
            a5,
            a6,
            a7,
            a8,
            a9,
            a10,
            a11,
            a12,
            a13,
            a14,
            a15,
            a16,
            a17,
            a18,
            a19,
            a20,
            a21,
            a22
          ),
          cvector = makeC(1.5, 13, 22),
          log = T,
          correction = 0
        ) +
          base::log(p.13),
        P_18 = ddirichletultinomial(
          nvector = base::c(
            chr1,
            chr2,
            chr3,
            chr4,
            chr5,
            chr6,
            chr7,
            chr8,
            chr9,
            chr10,
            chr11,
            chr12,
            chr13,
            chr14,
            chr15,
            chr16,
            chr17,
            chr18,
            chr19,
            chr20,
            chr21,
            chr22
          ),
          avector = base::c(
            a1,
            a2,
            a3,
            a4,
            a5,
            a6,
            a7,
            a8,
            a9,
            a10,
            a11,
            a12,
            a13,
            a14,
            a15,
            a16,
            a17,
            a18,
            a19,
            a20,
            a21,
            a22
          ),
          cvector = makeC(1.5, 18, 22),
          log = T,
          correction = 0
        ) +
          base::log(p.18),
        P_21 = ddirichletultinomial(
          nvector = base::c(
            chr1,
            chr2,
            chr3,
            chr4,
            chr5,
            chr6,
            chr7,
            chr8,
            chr9,
            chr10,
            chr11,
            chr12,
            chr13,
            chr14,
            chr15,
            chr16,
            chr17,
            chr18,
            chr19,
            chr20,
            chr21,
            chr22
          ),
          avector = base::c(
            a1,
            a2,
            a3,
            a4,
            a5,
            a6,
            a7,
            a8,
            a9,
            a10,
            a11,
            a12,
            a13,
            a14,
            a15,
            a16,
            a17,
            a18,
            a19,
            a20,
            a21,
            a22
          ),
          cvector = makeC(1.5, 21, 22),
          log = T,
          correction = 0
        ) +
          base::log(p.21)
      ) %>%
      dplyr::mutate(
        SumP = matrixStats::logSumExp(base::c(P_N, P_13, P_18, P_21)) %>%
          base::exp(),
        Posterior_N = base::exp(P_N) / SumP,
        Posterior_13 = base::exp(P_13) / SumP,
        Posterior_18 = base::exp(P_18) / SumP,
        Posterior_21 = base::exp(P_21) / SumP,
        maxP = base::max(base::c(
          Posterior_N,
          Posterior_13,
          Posterior_18,
          Posterior_21
        ))
      ) %>%
      dplyr::filter(SumP > 0) %>%
      dplyr::mutate(
        P_call = base::c(
          'No Aneuploidy',
          'Trisomy 13',
          'Trisomy 18',
          'Trisomy 21'
        )[base::which.max(base::c(P_N, P_13, P_18, P_21))]
      ) %>%
      dplyr::select(-dplyr::contains('alpha'))

    final_out_calls <- out_calls %>%
      dplyr::select(sample, total, P_call, maxP, dplyr::everything())
    base::
    return(dplyr::ungroup(final_out_calls))
  }
}

makeZscores <- function(
  indat,
  refType = 'auto',
  min_reads = 3e4,
  max_reads = 1e9
) {
  # Task: Calculates Z-scores per chromosome
  # Inputs:
  # - indat: full input tibble
  # - refType: "auto"="autosomal"; "sca"="sex chromosomal"; "diagnostic"
  # - min_reads: min number of reads per sample
  # - max_reads: max number of reads per sample
  # Checks made: yes

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

  base::
  return(dplyr::ungroup(z_out))
}

# IDs='EFE002'
# calls.auto=calls.auto.tf
# calls.sca=calls.sca.tf
# inDirichlet.auto <- inDirichlet.auto.tf
# inDirichlet.sca <- inDirichlet.sca.tf
# z.scores=z.scores.tf
# min_reads=4e4
# max_reads=1e9
plot.diagnostic <- function(
  IDs,
  calls.auto,
  calls.sca,
  z.scores,
  inDirichlet.auto,
  inDirichlet.sca,
  min_reads = 3e4,
  max_reads = 1e9,
  addLabels = FALSE
) {
  # Task: Plot three diagnostic plots
  # Inputs:
  # - IDs: a list or predefined piped string of sample IDs to look for
  # - calls.auto: the autosomal karyotype calls for the samples
  # - calls.sca: the sex chromosomal karyotype calls for the samples
  # - z.scores: the autosomal z-scores for the samples
  # - inDirichlet.auto: the dirichlet parameters for the autosomes
  # - inDirichlet.sca: the dirichlet parameters for the sex chromosomes
  # - min_reads: min number of reads per sample
  # - max_reads: max number of reads per sample
  # - addLabels: add "A", "B" and "C" to the plot panels?
  # Checks made: yes

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
    dplyr::filter(str_detect(sample, paste0(IDs, collapse = '|'))) %>%
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
      "minSamplesPerProtocol must be a positive integer, greater than ",
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
        dplyr::filter(str_detect(sample, paste0(IDs, collapse = '|'))),
      pch = 20,
      colour = 'black'
    ) +
    ggrepel::geom_label_repel(
      data = calls.sca %>%
        dplyr::select(sample, P_call, px, py) %>%
        dplyr::filter(str_detect(sample, paste0(IDs, collapse = '|'))),
      colour = 'black',
      box.padding = 1
    ) +
    ggsci::scale_colour_startrek(drop = FALSE, guide = 'none') +
    ggplot2::xlab('X-Rate') +
    ggplot2::ylab('Y-Rate') +
    ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x.bottom = ggplot2::element_blank(),
      axis.text.y.left = ggplot2::element_blank()
    )

  ## Second ggplot panel
  protocols <- z.scores %>%
    dplyr::filter(str_detect(sample, base::paste0(IDs, collapse = '|'))) %>%
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
  target_levels <- str_subset(
    z.scores$sample,
    base::paste0(IDs, collapse = '|')
  ) %>%
    base::unique() %>%
    base::unique() %>%
    base::c('Reference')
  gg2 <- z.scores %>%
    dplyr::filter(
      (Nj >= min_reads & Nj <= max_reads) |
        stringr::str_detect(sample, base::paste0(IDs, collapse = '|'))
    ) %>%
    dplyr::filter(protocol %in% protocols) %>%
    dplyr::mutate(
      isTarget = str_detect(sample, paste0(IDs, collapse = '|')),
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
        (str_detect(sample, base::paste0(IDs, collapse = '|')))
    ) %>%
    dplyr::filter(protocol %in% protocols) %>%
    dplyr::rename(sca_total = total) %>%
    dplyr::select(sample, protocol, P_call, sca_total, X, Y)
  gg3 <- dplyr::left_join(
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
        str_detect(chr, 'chr'),
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
        str_detect(chr, 'X|Y'),
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
    ggplot2::ggplot(ggplot2::aes(
      x = chr,
      y = fold_label,
      col = x_label,
      shape = x_label,
      label = base::round(fold)
    )) +
    ggplot2::theme_bw() +
    ggplot2::geom_hline(yintercept = 100, linetype = 'dashed') +
    ggplot2::geom_point(pch = 21, fill = 'black', size = 8) +
    ggplot2::geom_text(col = 'white', size = 3) +
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
    ggplot2::theme(legend.position = 'bottom')

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
  base::
  return(gg.final)
}

filterOutliers <- function(v) {
  # Task: filters for outliers
  # Hidden from user
  # Inputs:
  # - v: vector to filter

  if (!checkmate::checkDouble(v)) {
    base::stop('Input vector must be numeric.')
  }

  V <- EnvStats::rosnerTest(
    v,
    k = base::floor(base::length(v) * 0.05),
    warn = F
  )$all.stats %>%
    dplyr::filter(Outlier) %>%
    dplyr::pull(Obs.Num) %>%
    base::sort()
  v[V] <- NA
  base::
  return(v)
}

bound <- function(x) {
  # Task: bounds number to be in [0,1] (contamination estimation)
  # Hidden from user
  # Inputs:
  # - x: number
  if (!checkmate::checkDouble(x)) {
    base::stop('x must be a number.')
  }

  return(base::min(base::max(x, 0), 1))
}

adjustA <- function(a, c1, c2, g, correction = 1e-6) {
  # Adjusts alpha vector for contamination karyotypes
  # Hidden from user
  # Inputs:
  # - a: vector to change
  # - c1: X multiplier for karyotype 1
  # - c2: X multiplier for karyotype 2
  # - g: mixing parameter
  # - correction: error rate for Y mapping for non-Y karyotypes

  C1 <- base::c(c1, 1) %>%
    base::replace(. == 0, correction)
  C2 <- base::c(c2, 1) %>%
    base::replace(. == 0, correction)
  A1 <- (C1 * a) / base::sum(C1 * a) * base::sum(a)
  A2 <- (C2 * a) / base::sum(C2 * a) * base::sum(a)

  A <- (g * A1 + (1 - g) * A2) / base::sum(g * A1 + (1 - g) * A2)
  base::
  return(A)
}

estimateGamma <- function(nvector, avector, c1, c2, correction) {
  # Task: Estimate contamination from XX and XY based on observed counts
  # Hidden from user
  # Inputs:
  # - nvector: observed read counts
  # - avector: vector to change
  # - c1: X multiplier for karyotype 1
  # - c2: X multiplier for karyotype 2
  # - correction: error rate for Y mapping for non-Y karyotypes

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

  base::
  return(w1 * gamma1 + w2 * gamma2)
}

dirichlet.mle_filter <- function(y, whichK) {
  # Task: calculate Dirichlet parameters with filtering
  # Hidden from user
  # Inputs:
  # - y: input values
  # - whichK: which parameter to return

  x <- y[base::apply(y > 0, MARGIN = 1, all), ] %>%
    stats::na.omit()
  return(sirt::dirichlet.mle(x)$alpha[whichK])
}

ddirichletultinomial <- function(
  nvector,
  avector,
  cvector,
  log = T,
  correction = 1e-6
) {
  # Task: calculate Dirichlet probabilities with filtering
  # Hidden from user
  # Inputs:
  # - nvector: observed read counts
  # - avector: alpha vector to chnage
  # - cvector: change vector for new karyotype
  # - log: return value in log space
  # - correction: error rate for Y mapping for non-Y karyotypes

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
  base::
  return(dens)
}

makeC <- function(x, y, n) {
  # Task: make c vector for new karyotypes
  # Hidden from user
  # Inputs:
  # - x: value to replace 1
  # - y: location to make replacement
  # - n: length of vector

  X <- base::rep(1, n)
  if (y == 0) {
    base::
    return(X)
  } else {
    X[y] <- x
    base::
    return(X)
  }
}
