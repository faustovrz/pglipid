
#' @description The col files are tab delimited files representing info in data files
#' @details
#' @param input Biocyc col file
#'
#' @export
#' @return a dataframe. Note that many column names are serialized because
#'                      they are different values for a common key in dat files
#'
#' @examples
#' \dontrun{
#' read_col("proteins.dat")
#' }
#'
#'
read_col <- function(input = NULL) {
  read.table(
    input,
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    na.strings = "")
}

