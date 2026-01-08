utils::globalVariables(
  c(
    "Gamma",
    "P_13",
    "P_18",
    "P_21",
    "P_N",
    "P_X",
    "P_XX",
    "P_XXX",
    "P_XXY",
    "P_XY",
    "P_XYY",
    "P_call",
    "P_cont",
    "Posterior_13",
    "Posterior_18",
    "Posterior_21",
    "Posterior_N",
    "Posterior_X",
    "Posterior_XX",
    "Posterior_XXX",
    "Posterior_XXY",
    "Posterior_XY",
    "Posterior_XYY",
    "Posterior_cont",
    "SumC",
    "SumP",
    "a1",
    "a10",
    "a11",
    "a12",
    "a13",
    "a14",
    "a15",
    "a16",
    "a17",
    "a18",
    "a19",
    "a2",
    "a20",
    "a21",
    "a22",
    "a3",
    "a4",
    "a5",
    "a6",
    "a7",
    "a8",
    "a9",
    "correction",
    "maxP"
  )
)
#' Take prior and observed counts and calls karyotypes and Z-scores
#'
#' @param indat full input tibble
#' @param inDirichlet associated Dirichlet parameters
#' @param p_contamination the assumed probability of contamination
#'
#' @returns karyotype calls
#'
#' @export
callKaryotypes <- function(indat, inDirichlet, p_contamination = 0.1) {
  # Is it SCA or Auto?
  refType <- base::ifelse(
    any(stringr::str_detect(names(indat), '^px$')),
    'sca',
    'auto'
  )

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
    return(dplyr::ungroup(final_out_calls))
  }
}
