% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/networksDataPrep.R
\name{compileResultsMasterBL}
\alias{compileResultsMasterBL}
\title{LOAD BASELINE MASTER MATRICES (SPP X PIXEL MATRICES) FUNCTION}
\usage{
compileResultsMasterBL(
  file.ls,
  res.dir,
  out.dir,
  quant,
  dietcat,
  useCache = NULL,
  cacheRepo = options("reproducible.cachePath")
)
}
\arguments{
\item{file.ls}{is the list of file.paths for the master matrices}

\item{res.dir}{is the directory where scenario simulation results were saved}

\item{out.dir}{is the directory where outputs will be saved}

\item{quant}{is the quantile threshold of extinction chosen}

\item{dietcat}{diet categories (i.e., lowest trophic level nodes that are ubiquitous across all networks)}

\item{useCache}{is NULL, but will default to TRUE if the argument is not defined in the parent.frame()}

\item{cacheRepo}{passed to \code{reproducible::Cache}.}
}
\description{
compiles simulation results and saves the data.table
}
