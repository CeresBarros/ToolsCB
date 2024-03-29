% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/networksBetadiversity.R
\name{.networkTempBetaDiversityFamilies}
\alias{.networkTempBetaDiversityFamilies}
\title{NETWORK BETA DIVERSITY PER FAMILY USING WEBS AS INPUTS}
\usage{
.networkTempBetaDiversityFamilies(
  i,
  web1.ls,
  web2.ls,
  dietcat,
  charMatch,
  method,
  mode = "composition"
)
}
\arguments{
\item{i}{is the index to subset network lists by (e.g. pixel name if calculating temporal beta-diversity)}

\item{web1.ls}{is a list of networks to be indexed by i}

\item{web2.ls}{is a list of networks to be indexed by i}

\item{dietcat}{diet categories (i.e., lowest trophic level nodes that are ubiquitous across all networks)}

\item{charMatch}{is used to subset species in the networks beloging to a specific family.}

\item{method}{is passed to \code{networkTempBetaDiv}}

\item{mode}{is passed to \code{networkTempBetaDiv}}
}
\description{
Internal function used to calculate pairwise beta-diversity
for two networks per species family. The network pairs to compare should have the same
positions in the network lists supplied - e.g. the first network in list 1 is
compared with the first network in list 2, the second network in list 1 is
compared with the second network in list 2, and so on.
}
