#' CALCULATION OF THE MEAN FROM A BETA-INFLATED DISTRIBUTION
#'
#' adapted from `gamlss::meanBEINF`
#'
#' @param mu vector of values for the mu parameter Beta-inflated distribution
#' @param nu vector of values for the nu parameter Beta-inflated distribution
#' @param tau vector of values for the tau parameter Beta-inflated distribution
#'
#' @export
#'
#' @return expected mean Y

calcMeanBEINF <- function (mu, nu, tau) {
  meanofY <- (tau + mu)/(1 + nu + tau)
  meanofY
}
