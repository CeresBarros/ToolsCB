% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/networksDataPrep.R
\name{col_convert}
\alias{col_convert}
\title{CONVERT COLUMN NAMES IN MATRIX}
\usage{
col_convert(M, corresp, merge.fun = "max", na.rm = TRUE)
}
\arguments{
\item{M}{a \code{matrix} whose columns names will be changed}

\item{corresp}{is \code{data.frame} with rownames as the original column names
and a single column with the final column names. If more than one original
column matches a final one, the function will merge values using \code{merge.fun}}

\item{merge.fun}{funciton to use to merge columns. Default is \code{max},
which takes the maximum value across columns. Other possibilities are
\code{mean} or \code{min}.}

\item{na.rm}{is passed to merge.fun}
}
\description{
CONVERT COLUMN NAMES IN MATRIX
}
