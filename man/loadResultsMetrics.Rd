% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/networksDataPrep.R
\name{loadResultsMetrics}
\alias{loadResultsMetrics}
\title{LOAD METRICS RESULTS FUNCTION
compiles simulation resuls and saves the data.table}
\usage{
loadResultsMetrics(
  bl.dir = NULL,
  res.dir = NULL,
  out.dir = NULL,
  quant = c("No_ext_thresh", "min", "quant10", "quant25", "median", "quant50", "quant75",
    "quant90"),
  onlyBL = FALSE,
  onlyScen = FALSE
)
}
\arguments{
\item{bl.dir}{is the directory where baseline simulation results were saved}

\item{res.dir}{is the directory where scenario simulation results were saved}

\item{out.dir}{is the directory where outputs will be saved}

\item{quant}{is the quantile threshold of extinction chosen}

\item{onlyBL}{load only baseline networks?}

\item{onlyScen}{load only scenario networks?}
}
\description{
LOAD METRICS RESULTS FUNCTION
compiles simulation resuls and saves the data.table
}
