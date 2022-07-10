#' Modified `gamlss::getPEF`
#'
#' added `how = "mean"` and `level` to gamlss::getPEF
#' added `output` to chose between function or predicted values
#'
#' @param obj see `gamlss::getPEF`
#' @param term see `gamlss::getPEF`
#' @param data see `gamlss::getPEF`
#' @param n.points see `gamlss::getPEF`
#' @param parameter see `gamlss::getPEF`
#' @param level see `gamlss::getPEF`
#' @param type see `gamlss::getPEF`
#' @param how see `gamlss::getPEF`
#' @param fixed.at see `gamlss::getPEF`
#' @param plot see `gamlss::getPEF`
#' @param output chose between outputting the function or predicted values
#'
#' @export
#' @importFrom graphics abline
#' @importFrom stats splinefun median
#' @importFrom utils tail
getPEF.own <- function(obj = NULL, term = NULL, data = NULL, n.points = 100,
                       parameter = c("mu", "sigma", "nu", "tau"), level = NULL,
                       type = c("response", "link"), how = c("mean", "median", "last"), fixed.at = list(),
                       plot = FALSE, output = c("function", "vals")) {
  if (is.null(obj) || !class(obj)[1] == "gamlss")
    stop("Supply a standard GAMLSS model in obj")
  if (is.null(term))
    stop("The model term is not set")
  how <- match.arg(how)
  type <- match.arg(type)
  parameter <- match.arg(parameter)
  if (any(grepl("data", names(obj$call)))) {
    # DaTa <- if (startsWith(as.character(obj$call["data"]), "na.omit"))
    #   eval(parse(text = as.character(obj$call["data"])))
    # else get(as.character(obj$call["data"]))
    DaTa <- eval(parse(text = as.character(obj$call["data"])))
  } else if (is.null(data))
    stop("The data argument is needed in obj")

  DaTa <- as.data.frame(DaTa)

  mat <- matrix(0, nrow = dim(DaTa)[1] + n.points, ncol = dim(DaTa)[2])
  dat.temp <- as.data.frame(mat)
  names(dat.temp) <- v.names <- names(DaTa)
  pos <- which(names(dat.temp) == term)
  if (pos < 1)
    stop("supply a continuous term")
  if (is.factor(DaTa[, pos]))
    stop("the getPEF() is not suitable for factors")
  xvar <- seq(from = min(DaTa[, pos]), to = max(DaTa[, pos]),
              length.out = n.points)
  for (i in 1:dim(dat.temp)[2]) {
    if (pos == i) {
      dat.temp[, i] <- c(DaTa[, i], xvar)
    } else {
      ma <- fixed.at[[v.names[i]]]
      if (is.null(ma)) {
        if (how == "mean") {
          ma <- if (is.factor(DaTa[, i]))
            levels(DaTa[, i])[which.max(table(DaTa[, i]))]
          else mean(DaTa[, i])
        }
        if (how == "median") {
          ma <- if (is.factor(DaTa[, i]))
            levels(DaTa[, i])[which.max(table(DaTa[, i]))]
          else median(DaTa[, i])
        }
        if (how == "last") {
          ma <- if (is.factor(DaTa[, i]))
            levels(DaTa[, i])[which.max(table(DaTa[, i]))]
          else tail(DaTa[, i], 1)
        }
      }
      dat.temp[, i] <- c(DaTa[, i], rep(ma, n.points))
    }
  }
  fittted.orig <- predict(obj, newdata = tail(dat.temp, n.points),
                          type = type, parameter = parameter, level = level)
  theFun <- splinefun(xvar, fittted.orig)
  if (plot) {
    layout(matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE))
    plot(theFun(xvar) ~ xvar, ylab = "s()", xlab = term,
         type = "l")
    plot(theFun(xvar, deriv = 1) ~ xvar, xlab = term, ylab = "ds/dx",
         type = "l")
    abline(h = 0)
    layout(1)
  }
  invisible(theFun)
}
