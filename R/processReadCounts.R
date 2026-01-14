utils::globalVariables(
  c(
    "X",
    "Y",
    "auto",
    "chr1",
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
    "chr2",
    "chr20",
    "chr21",
    "chr22",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "n",
    "protocol",
    "total",
    "."
  )
)
#' processReadCounts
#'
#' @description
#' Task: takes read counts and makes useful file
#'
#'
#' @param rawReadCountsIn the reads counts via KP
#' @param refType for capture, use only on-target (F)?
#' @param minTotal minimum number of protocol counts to be included
#' @param minSamplesPerProtocol "auto"="autosomal"; "sca"="sex chromosomal"
#'
#' @returns outreadcounts
#'
#' @export
processReadCounts <- function(
  rawReadCountsIn,
  refType,
  minTotal = 1e3,
  minSamplesPerProtocol = 30
) {
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
      dplyr::select(tidyr::all_of(reqNames))
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
    dplyr::filter(dplyr::n() >= minSamplesPerProtocol) %>%
    dplyr::ungroup()

  # If refType is SCA then fold into X, Y and Automsomal
  if (tolower(refType) == "sca") {
    outReadCounts <- rawReadCounts %>%
      dplyr::mutate(
        auto = rowSums(dplyr::select(., tidyr::starts_with('chr')))
      ) %>%
      dplyr::select(sample, protocol, auto, X, Y) %>%
      dplyr::mutate(total = auto + X + Y) %>%
      dplyr::filter(total >= minTotal) %>%
      dplyr::mutate(px = X / total, py = Y / total, pz = auto / total)
  } else if (tolower(refType) == 'auto') {
    outReadCounts <- rawReadCounts %>%
      dplyr::select(-X, -Y) %>%
      dplyr::mutate(
        total = rowSums(dplyr::select(., tidyr::starts_with('chr')))
      ) %>%
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
  return(outReadCounts)
}
