% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatialTools.R
\name{neighboursMatrix}
\alias{neighboursMatrix}
\title{FIND NEIGHBOURS IN MATRIX}
\usage{
neighboursMatrix(mat)
}
\arguments{
\item{mat}{a \code{matrix}. Border of NA's is added to the matrix}
}
\value{
a matrix of 8 rows, with a columns per cell of the input matrix,
which is treated by columns
}
\description{
finds the 8 neighbours of each cell in a matrix
}
