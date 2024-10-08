% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/networkProperties.R
\name{calcSppTL}
\alias{calcSppTL}
\title{Calculate species trophic levels}
\usage{
calcSppTL(web = NULL, community = NULL, out.type = c("vector", "data.table"))
}
\arguments{
\item{web}{a square \code{matrix} representing an adjacency matrix.}

\item{community}{object of class "Community", output by \code{cheddar::Community}}

\item{out.type}{character. should trophic levels be output as a named vector (as the output of
\code{cheddar::PreyAveragedTrophicLevel}) or a data.table with species as columns
and a single row of trophic level values?}
}
\value{
a vector or \code{data.table} (see \code{out.type}) of trophic level values for each species in
the community\code{.}
}
\description{
Calculate species trophic levels
}
\details{
Cannibalistic links are removed with \verb{[cheddar::RemoveCannibalisticLinks()]}
before calculating species trophic levels with \verb{[cheddar::PreyAveragedTrophicLevel()]}
}
