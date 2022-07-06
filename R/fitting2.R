####################################################################
# Modified version of fitting from the brainwaver R package,
# removing unnecessary code beyond the fitting of a truncated power law
# and suppress histogram and writing outputs
####################################################################

#' @export
fitting2 <- function(degree.dist) {
  n.regions <- length(degree.dist)
  # tmp <- hist(degree.dist, breaks = c(0:nmax))
  # cum.dist <- 1 - cumsum(tmp$counts)/n.regions
  mu <- 1/(sum(degree.dist)/n.regions)
  nb <- length(degree.dist[degree.dist > 0])
  # gamma <- 1 + nb/(sum(log(degree.dist[degree.dist > 0])))
  x <- degree.dist
  x <- x[x > 0]
  n <- length(x)
  fn <- function(p) -(-n * p * log(sum(x)/(n * p)) - n * log(gamma(p)) +
                        (p - 1) * sum(log(x)) - n * p)
  out <- nlm(fn, p = 1, hessian = TRUE)
  alpha <- out$estimate
  beta <- sum(degree.dist)/(n.regions * alpha)
  #   AIC.exp <- -2 * (n.regions * log(mu) - mu * sum(degree.dist)) +
  #     2
  #   AIC.pow <- -2 * (n.regions * log(gamma - 1) - gamma * sum(log(x))) +
  #     2
  #   AIC.trunc <- -2 * (-out$minimum) + 2
  #   fitting <- "mu ="
  #   fitting <- paste(fitting, mu, sep = " ")
  #   fitting <- paste(fitting, "gamma = ", sep = "\n")
  #   fitting <- paste(fitting, gamma, sep = " ")
  #   fitting <- paste(fitting, "alpha = ", sep = "\n")
  #   fitting <- paste(fitting, alpha, sep = " ")
  #   fitting <- paste(fitting, "beta = ", sep = "\n")
  #   fitting <- paste(fitting, beta, sep = " ")
  #   fitting <- paste(fitting, "AIC exp = ", sep = "\n")
  #   fitting <- paste(fitting, AIC.exp, sep = " ")
  #   fitting <- paste(fitting, "AIC pow = ", sep = "\n")
  #   fitting <- paste(fitting, AIC.pow, sep = " ")
  #   fitting <- paste(fitting, "AIC trunc = ", sep = "\n")
  #   fitting <- paste(fitting, AIC.trunc, sep = " ")
  # write.table(fitting, "fitting.txt", row.names = FALSE, col.names = FALSE,
  # quote = FALSE)
  # list(mu = mu, gamma = gamma, alpha = alpha, beta = beta)
  list(alpha = alpha, beta = beta)
}
