#' Title
#'
#' @param rawReadCountsIn the reads counts for each chromosome
#' @param minSamplesPerProtocol minimum number of reads per protocol for parameter estimation
#' @param min_reads the minimum number of reads for Dirichlet parameter estimation
#' @param max_reads the maximum number of reads for Dirichlet parameter estimation
#' @param p_contamination the probability of a sample yielding significant contamination
#' @param show_plot show the clustering plot for sex chromosomal aneuploidy Dirichlet estimation?
#' @param printMissingIDs when combining karyotype calls, return names that are missing (or just the number of missing IDs)?
#'
#' @returns summary
#'
#' @export
#' @examples
#' runChASM(rawReadCountsIn = example_data)
#'
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
    # base::cat(
    #   "Warning: no protocol set, assuming all data to come from the same sequencing type (called 'default').\n"
    # )
    message(
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

  readcounts.auto <- processReadCounts(
    rawReadCountsIn = rawReadCountsIn,
    refType = 'auto',
    minSamplesPerProtocol = minSamplesPerProtocol
  )
  readcounts.sca <- processReadCounts(
    rawReadCountsIn = rawReadCountsIn,
    refType = 'sca',
    minSamplesPerProtocol = minSamplesPerProtocol
  )

  dirichlet.auto <- makeDirichlet(
    indat = readcounts.auto,
    refType = 'auto',
    min_reads = min_reads,
    max_reads = max_reads
  )
  dirichlet.sca <- makeDirichlet(
    indat = readcounts.sca,
    refType = 'sca',
    min_reads = min_reads,
    max_reads = max_reads,
    show_plot = show_plot
  )

  z.scores.auto <- makeZscores(
    indat = readcounts.auto,
    refType = 'auto',
    min_reads = min_reads,
    max_reads = max_reads
  )

  karyotypes.auto <- callKaryotypes(
    indat = readcounts.auto,
    inDirichlet = dirichlet.auto,
    p_contamination = p_contamination
  )
  karyotypes.sca <- callKaryotypes(
    indat = readcounts.sca,
    inDirichlet = dirichlet.sca,
    p_contamination = p_contamination
  )
  karyotypes <- combineData(
    calls.auto = karyotypes.auto,
    calls.sca = karyotypes.sca,
    z.scores = z.scores.auto,
    printMissingIDs = printMissingIDs
  )

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
