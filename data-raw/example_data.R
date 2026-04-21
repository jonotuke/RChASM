## code to prepare `example_data` dataset goes here
example_file <- fs::path_package("RChASM", "inst/extdata/example_data.tsv")
example_data <- readr::read_tsv(example_file)
usethis::use_data(example_data, overwrite = TRUE)
