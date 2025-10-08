## ----------------------------------------
## K-FOLD CROSS VALIDATION FUNCTIONS
## Ceres June 03 2020
## ----------------------------------------

#' CROSS-VALIDATION FUNCTION
#'
#' @param fullDT `data.table` with full dataset
#' @param statsModel the statistical model to validate. Only works with gamlss models
#' @param k integer with number of chunks that the data should be partitioned in
#' @param idCol column with pixel/observation IDs (optional)
#' @param sampleGroup column used as a grouping variable (i.e. random effect) in the  `statsModel`.
#'    Used to ensure that all folds contain all levels of this grouping variable.
#' @param origData the data used to fit the statsModule, needs to be passed to [`gamlss::predictAll()`]
#'   (it may not be able to access it) but also to make sure newdata in [`gamlss::predictAll()`]
#'  has the same variables (even if they're not used in the model)
#' @param level passed to `gamlss:::predict`
#' @param cacheObj1 object used by [`reproducible::Cache()`] for digesting,
#'  to avoid digesting the (potentially) large data arguments
#' @param cacheObj2  object used by [`reproducible::Cache()`] for digesting,
#'  to avoid digesting the (potentially) large data arguments
#' @param parallel logical. Uses [`future.apply::future_lapply()`] to parallelise
#'  model fitting across the k-folds, using `plan(multiprocess)`. Defaults to FALSE
#' @param cacheArgs a named `list` of arguments passed to inner `Cache` calls
#' @param classVar if a categorical (i.e. class) version of the response is available,
#'   it can be passed here to run validation metrics on class probabilities and accuracies.
#' @param ... further arguments passed to [`future::plan()`].
#'
#' @export
crossValidFunction <- function(fullDT, statsModel, origData, k = 4, idCol, sampleGroup = NULL,
                               parallel = FALSE, cacheArgs = NULL, level = NULL, classVar = NULL, ...) {

  dots <- list(...)

  if (!requireNamespace("gamlss", quietly = TRUE)) {
    stop("'gamlss' is not installed. Please install using:",
         "\ninstall.packages('gamlss')")
  }

  if (!is.null(idCol))
    origDataVars <- c(names(origData), idCol)

  ## remove NAs from the data without subsetting columns
  if (any(is.na(fullDT[, ..origDataVars])))
    stop("Please remove NAs from 'fullDT'")

  ## partition data into roughly equal chunks
  savedSeed <- .Random.seed
  on.exit(assign(".Random.seed", savedSeed, envir = .GlobalEnv), add = TRUE)
  set.seed(123)
  if (!is.null(sampleGroup)) {
    cols2 <- c(sampleGroup, idCol)
    sampDT <- fullDT[, ..cols2]
    sampDT[, sampID := sample(1:k, size = length(get(idCol)), replace = TRUE),
           by = get(sampleGroup)]
    rm(cols2)
  } else {
    sampDT <- unique(fullDT[, ..idCol])
    sampDT[, sampID := sample(1:k, size = length(get(idCol)), replace = TRUE)]
  }

  ## join samp IDs with data
  fullDT <- sampDT[fullDT, on = idCol]

  origDataVars <- c(origDataVars, "sampID")

  message(paste("Starting cross-validation using", k, "folds"))

  if (parallel) {
    if (!requireNamespace("doSNOW", quietly = TRUE)) {
      stop("'doSNOW' is not installed. Please install using:",
           "\ninstall.packages('doSNOW')")
    }

    if (!requireNamespace("foreach", quietly = TRUE)) {
      stop("'foreach' is not installed. Please install using:",
           "\ninstall.packages('foreach')")
    } else require(foreach)

    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("'parallel' is not installed. Please install using:",
           "\ninstall.packages('parallel')")
    }

    # Extract relevant arguments
    clusterArgs <- formalArgs(makeCluster)
    clusterArgs <- dots[intersect(names(dots), clusterArgs)]

    # Set up cluster based on OS
    if (Sys.info()[["sysname"]] == "Windows") {
      cl <- do.call(makeCluster, clusterArgs)
      doSNOW::registerDoSNOW(cl)
    } else {
      clusterArgs$type <- "FORK"
      cl <- do.call(makeCluster, clusterArgs)
      doSNOW::registerDoSNOW(cl)
    }
    on.exit(stopCluster(cl), add = TRUE)

    # Set up progress bar
    pb <- txtProgressBar(min = 0, max = length(unique(fullDT$sampID)[1:2]), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    pkgs <- c("caret",
              "data.table",
              "gamlss",
              "gamlss.inf",
              "reproducible",
              "stats")

    # Apply function in parallel using foreach, move some arguments to a separate lise
    crossValidResults <- foreach(sampID = unique(fullDT$sampID)[1],
                                 .packages = pkgs,
                                 .options.snow = opts) %dopar% {
                                   calcCrossValidMetrics(samp = sampID,
                                                         classVar = classVar,
                                                         idCol = idCol,
                                                         fullDT = fullDT,
                                                         origData = origData,
                                                         statsModel = statsModel,
                                                         origDataVars = origDataVars,
                                                         level = level,
                                                         cacheArgs = cacheArgs)
                                 }

  } else {
    crossValidResults <- unique(fullDT$sampID) |>
      purrr::map(\(x) calcCrossValidMetrics(x,
                                            classVar = classVar,
                                            idCol = idCol,
                                            fullDT = fullDT,
                                            origData = origData,
                                            statsModel = statsModel,
                                            origDataVars = origDataVars,
                                            level = level,
                                            cacheArgs = cacheArgs))
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
#' @param classVar if a categorical (i.e. class) version of the response is available,
#'   it can be passed here to run validation metrics on class probabilities and accuracies.
#' @param idCol row IDs. Needed when `!is.null(classVar)`
#' @param origDataVars a character vector of the variables used in model fitting (including response variable and random effects.)
#' @param cacheArgs a named `list` of arguments passed to `Cache`
#'
#' @return a list with 2 entries
#'
#' @importFrom reproducible Cache
#' @importFrom data.table as.data.table set
#' @importFrom stats update
#' @importFrom gamlss predict Rsq
#'
#' @export
calcCrossValidMetrics <- function(samp, fullDT, origData, statsModel,
                                  origDataVars, level = NULL, cacheArgs = NULL,
                                  classVar = NULL, idCol = NULL) {
  if (!requireNamespace("gamlss", quietly = TRUE)) {
    stop("'gamlss' is not installed. Please install using:",
         "\ninstall.packages('gamlss')")
  }

  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("'caret' is not installed. Please install using:",
         "\ninstall.packages('caret')")
  }

  message(paste("Fold", samp))
  ## predict requires the original and new data to have the same columns
  if (!all(names(origData) %in% names(fullDT)))
    stop("'fullDT' needs to include all the columns in 'origData'")

  ## subset
  trainData <<- fullDT[sampID != samp, ..origDataVars]
  testData <- fullDT[sampID == samp, ..origDataVars]

  ## checks, if there's a RE variable, all it's levels need to exist in the sample
  REvar <- names(statsModel$mu.coefSmo[[1]]$coefficients$random)
  if (!is.null(REvar)) {
    if (!REvar %in% colnames(fullDT)) {
      stop("Grouping variable not found in 'fullDT'")
    }

    if (length(setdiff(unique(fullDT[[REvar]]),
                       unique(testData[[REvar]]))) |
        length(setdiff(unique(fullDT[[REvar]]),
                       unique(trainData[[REvar]]))))
      stop("Fires lost in sampling!")
  }


  if (any(is.na(trainData)) | any(is.na(testData)))
    stop("Please remove NAs from the variables going in the model")

  ## trainData an testData cannot have extra cols with respect to those in the original
  ## data used to fit the model
  cols <- names(origData)
  trainData <<- trainData[, ..cols]
  testData <- testData[, ..cols]

  ## test
  # trainData <<- trainData[1:200, ..cols]
  # testData <- testData[201:301, ..cols]

  cacheObj <- list(resid(statsModel), samp)

  ## refit model on training sample then predict
  trainModel <- tryCatch({
    update(object = statsModel,
            data = trainData) #|>
  # Cache(omitArgs = c("data", "object"),   ## 1) need to pipe this cache call, gamlss won't work, 2) omitArgs is not working properly here, but all seems good in global env.
        # .cacheExtra = cacheObj)
  }, error = function(e) e)

  if (is(trainModel, "error")) {
    message("Model could not be re-fit. Error:")
    message(trainModel)
    validMetrics <- c("RMSE" = NA, "Rsquared" = NA, "MAE" = NA,
                      "Rsq" = NA, "TGD" = NA, "predictError" = NA)
    return(validMetrics)
  }

  params <- c("mu", "nu", "tau")
  predictionsDT <- sapply(params, FUN = function(param) {
    predict(object = trainModel,
            what = param,
            newdata = testData,
            data = trainData,
            type = "response",
            level = level) |>
      Cache(omitArgs = c("object", "newdata", "data"),
            .cacheExtra = c(cacheObj, param))
  }) |>
    as.data.table()

  ## add response variable
  modform <- formula(statsModel)
  modterms <- terms(modform)
  respVar <- as.character(attr(modterms, "variables")[attr(modterms, "response") + 1])

  set(predictionsDT, NULL, "obs", testData[, get(respVar)])

  ## predict using meanBEINF approach
  if (trainModel$family[1] != "BEINF")
    stop("the object does not have a BEINF distribution")

  predictionsDT[, pred := calcMeanBEINF(mu, nu, tau)]

  ## VALIDATION STATISTICS WITH CLASSES -----------------------
  if (!idCol %in% colnames(testData)) {
    testData <- na.omit(fullDT[sampID == samp, ..origDataVars]) ## redo testData in case idCol was dropped when subsetting to model data
  }

  if (exists("classVar")) {
    ## add severity classes
    predictionsDT[, c(idCol) := testData[[idCol]]]
    cols <- c(idCol, classVar)
    predictionsDT <- fullDT[, ..cols][predictionsDT, on = idCol]

    ## convert to classes, using the quantiles corresponding to the observed class proportions
    ## accumulate proportions to get probabilities
    quantProbs <- cumsum(table(predictionsDT[[classVar]])/nrow(predictionsDT))
    classRanges <- c(0, quantile(predictionsDT$pred, probs = quantProbs))

    predictionsDT[, predCLASS := cut(pred, breaks = classRanges,
                                     include.lowest = TRUE, right = FALSE)]  ## classify as with intervals as ],]

    ## convert to numbered factor (subtracting one, because classes are 0-5)
    predictionsDT[, predCLASS := as.numeric(predCLASS)-1]
    classes <- as.character(sort(unique(fullDT[[classVar]])))
    predictionsDT[, `:=`(obsCLASS = factor(get(classVar), levels = classes),
                         predCLASS = factor(predCLASS, levels = classes))]

    ## VALIDATION STATISTICS WITH CLASSES ----------------------------------
    ## calculate overall statistics
    validMetricsClass <- caret::multiClassSummary(predictionsDT[, list(obs = obsCLASS,
                                                                       pred = predCLASS)],
                                                  lev = classes)
    ## calculate confusion matrix
    confMatrix <- caret::confusionMatrix(data = predictionsDT$predCLASS,
                                         reference = predictionsDT$obsCLASS)

    out <- list(validMetricsClass = validMetricsClass,
                confMatrix = confMatrix)
  }

  ## VALIDATION STATISTICS WITH CONTINUOUS VARIABLE -----------------------
  RsqGAMLSS <- Rsq(object = trainModel) |>
    Cache(omitArgs = c("object"),
          .cacheExtra = cacheObj)
  TGDstats <- getTGD(object = trainModel,
                     newdata = testData,
                     data = trainData) |>
    Cache(omitArgs = c("object", "newdata", "data"),
          .cacheExtra = cacheObj)

  caretSumm <- caret::defaultSummary(data.frame(obs = predictionsDT$obs, pred = predictionsDT$pred))
  RMSE <- caret::RMSE(pred = predictionsDT$pred, obs = predictionsDT$obs)
  validMetricsCont <- c(caretSumm,
                        RsqGAMLSS = RsqGAMLSS,
                        RMSE = RMSE,
                        TGD = TGDstats$TGD,
                        predictError = TGDstats$predictError)

  if (!exists("out"))
    out <- list()

  out <- append(out,
                list(
                  validMetrics = validMetricsCont,
                  validMetricsClass = validMetricsClass,
                  confMatrix = confMatrix,
                  coefs = gamlss::coefAll(trainModel)
                )
  )

  return(out)
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

#' Root Mean Squared Error
#'
#' (Modified from `Metrics::rmse`)
#'
#' \code{rmse} computes the root mean squared error between two numeric vectors
#'
#' @param actual The ground truth numeric vector.
#' @param predicted The predicted numeric vector, where each element in the vector
#'                  is a prediction for the corresponding element in \code{actual}.
#' @author Michael Frasco
#' @examples
#' actual <- c(1.1, 1.9, 3.0, 4.4, 5.0, 5.6)
#' predicted <- c(0.9, 1.8, 2.5, 4.5, 5.0, 6.2)
#' rmse(actual, predicted)
rmse <- function(actual, predicted) {
  return(sqrt(mean((actual - predicted) ^ 2)))
}
