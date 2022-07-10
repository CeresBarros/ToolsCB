## ----------------------------------------
## K\-FOLD CROSS VALIDATION FUNCTIONS
## Ceres June 03 2020
## ----------------------------------------

#' CROSS\-VALIDATION FUNCTION
#'
#' @param fullDT `data.table` with full dataset
#' @param statsModel the statistical model to validate. Only works with gamlss models
#' @param k integer with number of chunks that the data should be partitioned in
#' @param idCol column with pixel/observation IDs (optional)
#' @param origData the data used to fit the statsModule, needs to be passed to `gamlss::predictAll`
#'   (it may not be able to access it) but also to make sure newdata in `gamlss::predictAll`
#'  has the same variables (even if they're not used in the model)
#' @param level passed to `gamlss:::predict`
#' @param cacheObj1 object used by `reproducible::Cache` for digesting,
#'  to avoid digesting the (potentially) large data arguments
#' @param cacheObj2  object used by `reproducible::Cache` for digesting,
#'  to avoid digesting the (potentially) large data arguments
#' @param parallel logical. Uses `future.apply::future_lapply` to parallelise
#'  model fitting across the k-folds, using `plan(multiprocess)`. Defaults to FALSE
#' @param ... further arguments passed to `future::plan`
#' @param cacheArgs a named `list` of arguments passed to inner `Cache` calls
#'
#' @importFrom future.apply future_lapply
#' @importFrom future plan
#' @export
crossValidFunction <- function(fullDT, statsModel, origData, k = 4, idCol,
                                parallel = FALSE, cacheObj1 = NULL, cacheObj2 = NULL,
                                cacheArgs = NULL, level = NULL, ...) {
  if (!is.null(idCol))
    origDataVars <- c(names(origData), idCol)

  ## remove NAs from the data without subsetting columns
  if (any(is.na(fullDT[, ..origDataVars])))
    stop("Please remove NAs from the variables going in the model")

  ## partition data into roughly equal chunks
  sampDT <- unique(fullDT[, ..idCol])
  sampDT[, sampID := sample(1:k, size = length(get(idCol)), replace = TRUE)]
  ## join samp IDs with data
  fullDT <- sampDT[fullDT, on = idCol]

  origDataVars <- c(origDataVars, "sampID")

  message(paste("Starting cross-validation using", k, "folds"))
  if (parallel) {
    if (Sys.info()[["sysname"]] == "Windows") {
      plan(multisession, gc = TRUE, ...)
    } else plan(multicore, ...)
    crossValidResults <- future_lapply(unique(fullDT$sampID), FUN = calcCrossValidMetrics,
                                       fullDT = fullDT, origData = origData,
                                       statsModel = statsModel, origDataVars = origDataVars,
                                       level = level, cacheArgs = cacheArgs)
    ## Explicitly close workers
    future:::ClusterRegistry("stop")
  } else {
    crossValidResults <- lapply(unique(fullDT$sampID), FUN = calcCrossValidMetrics,
                                fullDT = fullDT, origData = origData,
                                statsModel = statsModel, origDataVars = origDataVars,
                                level = level, cacheArgs = cacheArgs)
  }
  return(crossValidResults)
}


#' CALCULATE VALIDATION METRICS AND CONFUSION MATRIX
#'
#' to allow caching without digesting the large data table
#'
#' @param samp the sample number to pick to use as the test data set
#' @param fullDT the full dataset (not necessarily the one used to fit
#'   `statsModel`, which could have been a subset (e.g. fewer columns))
#' @param origData the data used to fit `statsModel`
#' @param statsModel the fitted model
#' @param level passed to `gamlss:::predict`
#' @param origDataVars a character vector of the variables used in model fitting (including response variable and random effects.)
#' @param cacheArgs a named `list` of arguments passed to `Cache`
#'
#' @return a list with 2 entries
#'
#' @importFrom gamlss Rsq getTGD
#' @importFrom reproducible Cache
#' @importFrom data.table as.data.table set
#' @importFrom caret defaultSummary
#' @importFrom stats update
#'
#' @export
calcCrossValidMetrics <- function(samp, fullDT, origData, statsModel, origDataVars, level = NULL, cacheArgs = NULL) {
  message(paste("Fold", samp))
  ## predict requires the original and new data to have the same columns
  if (!all(names(origData) %in% names(fullDT)))
    stop("'fullDT' needs to include all the columns in 'origData'")

  ## subset
  trainData <<- fullDT[sampID != samp, ..origDataVars]
  testData <- fullDT[sampID == samp, ..origDataVars]

  if (any(is.na(trainData)) | any(is.na(testData)))
    stop("Please remove NAs from the variables going in the model")

  ## trainData an testData cannot have extra cols with respect to those in the original
  ## data used to fit the model
  cols <- names(origData)
  trainData <<- trainData[, .SD, .SDcols = cols]   ## need to export to .Global for gamlss...
  testData <- testData[, .SD, .SDcols = cols]

  ## refit model on training sample then predict
  trainModel <- tryCatch(update(object = statsModel, data = trainData), error = function(e) e)

  if (is(trainModel, "error")) {
    message("Model could not be re-fit. Error:")
    message(trainModel)
    validMetrics <- c("RMSE" = NA, "Rsquared" = NA, "MAE" = NA,
                      "Rsq" = NA, "TGD" = NA, "predictError" = NA)
  } else {
    params <- c("mu", "nu", "tau")
    names(params) <- params
    predictionsDT <- lapply(params, FUN = function(param) {
      predict(trainModel, what = param,
              newdata = testData, data = trainData,
              type = "response", level = level)
    })
    predictionsDT <- as.data.table(do.call(cbind, predictionsDT))

    ## add response variable
    set(predictionsDT, NULL, "invRobust", testData$invRobust)

    ## predict using meanBEINF approach
    if (trainModel$family[1] != "BEINF")
      stop("the object does not have a BEINF distribution")

    predictionsDT[, predinvRobust := .calcMeanBEINF(mu, nu, tau),
                  by = row.names(predictionsDT)]

    ## VALIDATION STATISTICS WITH CONTINUOUS VARIABLE -----------------------
    RsqGAMLSS <- Rsq(trainModel)
    TGDstats <- getTGD(trainModel, newdata = testData, data = trainData)

    validMetrics <- c(defaultSummary(data.frame(obs = predictionsDT$invRobust, pred = predictionsDT$predinvRobust)),
                      "Rsq" = RsqGAMLSS,
                      TGD = TGDstats$TGD,
                      predictError = TGDstats$predictError)
  }
  list(validMetrics = validMetrics)
}

#' CACHE-COMPATIBLE MODEL UPDATE
#'
#' to allow caching without digesting the large data table and model
#'
#' @param ... arguments passed to `update`
#' @param cacheObj1 an object used by Cache for digesting, to avoid digesting
#'   the (potentially) large data arguments e.g.: model coefficients and a
#'   sample of a column drawn with a set seed.
#' @param cacheObj2 an object used by Cache for digesting, to avoid digesting
#'   the (potentially) large data arguments e.g.: model coefficients and a
#'   sample of a column drawn with a set seed.
#'
#' @export
updateModelCached <- function(..., cacheObj1 = NULL, cacheObj2 = NULL) {
  updatedModel <- update(...)
  return(updatedModel)
}


#' CALCULATION OF THE MEAN FROM A BETA-INFLATED DISTRIBUTION
#' adapted from `gamlss::meanBEINF`
#' @param mu vector of values for the mu parameter of
#'   the Beta-inflated distribution
#' @param nu vector of values for the nu parameter of
#'   the Beta-inflated distribution
#' @param tau vector of values for the tau parameter of
#'   the Beta-inflated distribution
.calcMeanBEINF <- function(mu, nu, tau) {
  meanofY <- (tau + mu)/(1 + nu + tau)
  meanofY
}
