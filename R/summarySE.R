#' Summarize data
#'
#' Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95\\%).
#'
#' @param data a data frame.
#' @param measurevar the name of a column that contains the variable to be summarized
#' @param groupvars a vector containing names of columns that contain grouping variables
#' @param na.rm a boolean that indicates whether to ignore NA's, defaults to FALSE.
#' @param conf.interval the percent range of the confidence interval (default is 95\\%)
#' @param .drop passed to `ddply`
#'
#' @export
#'
#' @return A `data.frame` of summary statistics of each variable
#'  by grouping variables
#'
#' @importFrom plyr ddply rename
#' @importFrom stats qt
#' @importFrom stats sd var

summarySE <- function(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE,
                      conf.interval = 0.95, .drop = TRUE) {

  ## This does the summary. For each group's data frame, return a vector with
  ## N, mean, and sd
  datac <- ddply(data, groupvars, .drop = .drop,
                 .fun = function(xx, col) {
                   c(N = length2(xx[[col]], na.rm = na.rm),
                     mean = mean(xx[[col]], na.rm = na.rm),
                     sd = sd(xx[[col]], na.rm = na.rm),
                     var = var(xx[[col]], na.rm = na.rm))
                 },
                 measurevar)

  ## Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))

  datac$se <- datac$sd / sqrt(datac$N)  ## Calculate standard error of the mean

  ## Confidence interval multiplier for standard error
  ## Calculate t-statistic for confidence interval:
  ## e.g., if conf.interval is 0.95, use 0.975 (above/below), and use df = N-1
  ciMult <- qt(conf.interval/2 + 0.5, datac$N - 1)
  datac$ci <- datac$se * ciMult

  return(datac)
}


#' New version of length which can handle NA's
#'
#' @param x a vector
#' @param na.rm if TRUE does not acount for NA's when taking the length of `x`.
#'   Defaults to FALSE which has the same behaviour as `lenght`.
#'
#' @export
#'
#' @return an integer of the length of `x` accounting for NAs or not.

length2 <- function(x, na.rm = FALSE) {
  if (na.rm) sum(!is.na(x)) else length(x)
}
