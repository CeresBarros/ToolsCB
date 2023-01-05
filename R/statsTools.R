#' CALCULATION OF THE MEAN FROM A ZERO- AND ONE-INFLATED BETA-INFLATED DISTRIBUTION
#'
#' adapted from `gamlss::meanBEINF`
#'
#' @param mu vector of values for the mu parameter of the zero- and one-inflated Beta distribution (BEINF)
#' @param nu vector of values for the nu parameter of the BEINF
#' @param tau vector of values for the tau parameter of the BEINF
#'
#' @export
#'
#' @return expected mean Y

calcMeanBEINF <- function (mu, nu, tau) {
  meanofY <- (tau + mu)/(1 + nu + tau)
  return(meanofY)
}

#' CALCULATION OF THE PROBABILITY OF ZERO FROM A ZERO- AND ONE-INFLATED BETA-INFLATED DISTRIBUTION
#'
#' Follows the parameterisation used in `gamlss` package. 
#' See Rigby et al. 2020 (Distributions for Modeling Location, Scale, and Shape: using GAMLSS in R)
#'
#' @param nu vector of values for the nu parameter of the zero- and one-inflated Beta distribution (BEINF)
#' @param tau vector of values for the tau parameter of the BEINF
#'
#' @export
#'
#' @return expected probability of 0 [P(Y=0)]
calcP0BEINF(nu, tau) {
  p0 <- nu/(1 + nu + tau)
  return(p0)
}

#' CALCULATION OF THE PROBABILITY OF ONE FROM A ZERO- AND ONE-INFLATED BETA-INFLATED DISTRIBUTION
#'
#' Follows the parameterisation used in `gamlss` package. 
#' See Rigby et al. 2020 (Distributions for Modeling Location, Scale, and Shape: using GAMLSS in R)
#'
#' @param nu vector of values for the nu parameter of the zero- and one-inflated Beta distribution (BEINF)
#' @param tau vector of values for the tau parameter of the BEINF
#'
#' @export
#'
#' @return expected probability of 1 [P(Y=1)]
calcP1BEINF(nu, tau) {
  p1 <- tau/(1 + nu + tau)
  return(p1)
}

#' CALCULATION OF THE VARIANCE OF A ZERO- AND ONE-INFLATED BETA-INFLATED DISTRIBUTION
#'
#' Follows the parameterisation used in `gamlss` package, using variance estimation described in 
#'  and Ferrari (2010). See also Rigby et al. 2020 (Distributions for Modeling Location,
#'  Scale, and Shape: using GAMLSS in R)
#' 
#' @param mu vector of values for the mu parameter of the zero- and one-inflated Beta distribution (BEINF)
#' @param sigma vector of values for the sigma parameter of the BEINF
#' @param nu vector of values for the nu parameter of the BEINF
#' @param tau vector of values for the tau parameter of the BEINF
#'
#' @export
#'
#' @return expected variance of Y
calcVarBEINF(mu, sigma, nu, tau) {
  gamma <- (tau*(1 + nu + tau))/((nu + tau)*(1 + nu + tau))
  alpha <- (nu + tau)/(1 + nu + tau)
  Vmu <- mu*(1 - mu)
  V1 <- gamma*(1 - gamma)
  V2 <- Vmu/(sigma + 1)

  VarY <- alpha*V1 + (1 - alpha)*V2 + alpha*(1 - alpha)(gamma - mu)^2
  return(VarY)
}