#' A function for saving the results of a ChASM analysis as a tsv
#'
#' @param inChASM the result of a full ChASM analysis (output from runChASM)
#' @param file the full path and file name to place output
#' @param sort_by_samplename reorder alphabetically by sample names?
#'
#' @returns saves results
#'
#' @export
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
