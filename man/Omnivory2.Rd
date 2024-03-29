% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_nu17May_Ceres.R
\name{Omnivory2}
\alias{Omnivory2}
\title{Network omnivory}
\usage{
Omnivory2(M, dietcat = NULL, level = PreyAveragedTrophicLevel)
}
\arguments{
\item{M}{a square \code{matrix} representing an adjacency newtork}

\item{dietcat}{diet categories (i.e., lowest trophic level nodes that are ubiquitous across all networks)}

\item{level}{a function to calculate the trophic level. See \code{cheddar}
for options}
}
\description{
Network omnivory
}
