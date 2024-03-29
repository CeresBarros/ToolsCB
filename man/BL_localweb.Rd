% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/networksBuildLocalWebs.R
\name{BL_localweb}
\alias{BL_localweb}
\title{Function to calculate local webs - Baseline}
\usage{
BL_localweb(
  metaweb = NULL,
  SPPCODE = NULL,
  SPP.HAB = NULL,
  PIX.HAB = NULL,
  dietcat = NULL,
  HELP = TRUE
)
}
\arguments{
\item{metaweb}{a square \code{matrix} representing an adjacency metaweb}

\item{SPPCODE}{one or more species IDs}

\item{SPP.HAB}{a spp x habitat \code{matrix}}

\item{PIX.HAB}{a pix x habitat \code{matrix}}

\item{dietcat}{diet categories (i.e., lowest trophic level nodes that are ubiquitous across all networks)}

\item{HELP}{print help}
}
\description{
Function to calculate local webs - Baseline
}
