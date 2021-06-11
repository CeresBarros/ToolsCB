#' REPLACE NAs BY ANOTHER VALUE IN A DATA.TABLE
#'
#' @param x is a \code{data.table}
#' @param val is the value to replace NAs for. By default, it looks
#'  for and replaces 0s
#'
#' @export
#'
#' @return a \code{data.table}
#'
#' @importFrom data.table data.table

replaceNAs <- function(x, val = 0) {
  x[is.na(x)] <- val
  x
}
