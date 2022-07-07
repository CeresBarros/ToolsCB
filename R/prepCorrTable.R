#' LATEX/RMARKDOWN CORRELATION TABLES
#' @description  function to make a nice correlation table for latex/Rmd
#'               usage example:
#'               pTable <- rcorr()
#'               tabl <- prepCorrTable(DT) %>% xtable::xtable(.)
#'               print(tabl, type = "html.
#'               Adapted from http://myowelt.blogspot.com/2008/04/beautiful-correlation-tables-in-r.html
#' @param x a matrix or any object compatible with as.matrix, from where the covariance matrix will be calculated.
#'          alternatively, corTable and pTable matrices can be supplied.
#' @param method see \code{Hmisc::rcorr} \code{type} argument.
#' @param corTable matrix of correlation values obtained from e.g. \code{Hmisc::rcorr(...)$r}.
#'                 If supplied, x will be ignored
#' @param pTable matrix of pTable-values values obtained from e.g. \code{Hmisc::rcorr(...)$pTable}
#'
#' @export
#' @importFrom Hmisc rcorr

prepCorrTable <- function(x = NULL, corTable = NULL, pTable = NULL, method = "pearson") {
  ## checks
  if (!is.null(x)) {
    if (!all(class(x) == "matrix"))
      x <- as.matrix(x)
  } else if (is.null(corTable) | is.null(pTable)) {
    stop("x not supplied, please provide 'corTable' and 'pTable'")

    if (is.null(corTable) != is.null(pTable)) {
      stop("Please provide both 'corTable' and 'pTable'")
    }
  }

  if (!method %in% c("pearson", "spearman"))
    stop("method must be 'pearson' or 'spearman'")

  if (is.null(corTable) & is.null(pTable)) {
    corTable <- rcorr(x, type = method)$r
    pTable <- rcorr(x, type = method)$P
  }

  ## define notions for significance levels; spacing is important.
  mystars <- ifelse(pTable < .001, "***", ifelse(pTable < .01, "** ", ifelse(pTable < .05, "* ", " ")))

  ## trunctuate the matrix that holds the correlations to two decimal
  corTable <- format(round(cbind(rep(-1.11, ncol(corTable)), corTable), 2))[,-1]

  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(corTable, mystars, sep=""), ncol=ncol(corTable))
  diag(Rnew) <- paste(diag(corTable), " ", sep="")
  rownames(Rnew) <- colnames(corTable)
  colnames(Rnew) <- paste(colnames(corTable), "", sep="")

  ## remove upper triangle
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew)

  ## remove last column and return the matrix (which is now a data frame)
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  return(Rnew)
}
