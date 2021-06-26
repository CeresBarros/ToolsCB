#' CALCULATE AND COMPARE TWO HYPERVOLUMES
#'
#' Calculates two \code{hypervolumes} on two sets of "raw" data or on the factor scores of
#' an ordination done on the entire dataset (across the two sets of data).
#'
#' @param HVdata1 the \code{data.frame} from which the first hypervolume will be calculated
#' @param HVdata2 the \code{data.frame} from  which the second hypervolume will be calculated.
#'   \code{HVdata1} and \code{HVdata2} must have the same column names and ideally the same number of rows.
#' @param ordination determines whether hypervolumes will be calculated on the raw data ("none") of on
#'   ordination factor scores (one of "PCA", "HillSmith", "dudi.mix"). Defaults to using a PCA.
#' @param HVidvar the name or number of the column (in \code{HVdata1/2}) containing the ID of each hypervolume.
#' @param init.vars is a vector with the variable indices (column numbers) to use. Defaults to NULL
#'   and all variables are used - ATTENTION the user should not have more than 8 variables unless
#'   using an ordination to reduce dimensionality. See Blonder et al 2014.
#' @param noAxes determines the number of axes/columns to use for hypervolume calculation
#' @param do.scale activates scaling prior to the ordination/calculation of HVs
#' @param freeBW determines whether a bandwidth estimator will be used to calculate bandwith per variable.
#'    Only used if \code{HVmethod} is "box" or "gaussian".
#' @param bwHV1 determines the bandwidth value for first hypervolume. Only used if \code{freeBW} is FALSE. see \code{hypervolume::hypervolume}
#'    Only used if \code{HVmethod} is "box" or "gaussian".
#' @param bwHV2 determines the bandwidth value for second hypervolume. Only used if \code{freeBW} is FALSE. see \code{hypervolume::hypervolume}
#'    Only used if \code{HVmethod} is "box" or "gaussian".
#' @param HVmethod determines the method used to calculate hypervolumes - passed to \code{hypervolume::hypervolume} method argument
#' @param no.runs determines how many times the HV calculations and comparisons are repeated
#' @param plotOrdi activates/desactivates plotting for ordinations - plots are saved to PDFs not plotted interactively
#' @param plotHV activates/desactivates plotting for hypervolumes - plots are saved to PDFs not plotted interactively
#' @param saveOrdi activates/desactivates saving ordinations outputs
#' @param outputs.dir is the directory to store results
#' @param file.suffix is a character string used as a suffix in file names
#' @param verbose is passed to \code{hypervolume::*} functions to control printing of diagnostic outputs
#' @param ... further arguments passed to \code{hypervolume::hypervolume}.
#'
#' @return Nothing is returned. Hypervolumes and comparisons are saved \code{outputs.dir}
#'
#' @export
#'
#' @import data.table
#' @import pdist
#' @importFrom grDevices cairo_pdf graphics.off
#' @importFrom utils write.table

hypervolumes <- function(HVdata1, HVdata2, HVidvar, ordination = "PCA", init.vars = NULL,
                         noAxes = NULL, do.scale = FALSE, HVmethod = "box",
                         freeBW = FALSE, bwHV1 = NULL, bwHV2 = NULL,
                         no.runs = 1, plotOrdi = TRUE, plotHV = TRUE, saveOrdi = TRUE,
                         outputs.dir, file.suffix, verbose = TRUE, ...) {
  ## do some checks:
  if (!all(is(HVdata1, "data.frame"), is(HVdata2, "data.frame")))
    stop("HVdata1 and HVdata2 are not two dataframes")
  if (length(setdiff(names(HVdata1), names(HVdata2))) |
      length(setdiff(names(HVdata2), names(HVdata1))))
    stop("HVdata1 and HVdata2 columns do not match")
  if (!ordination %in% c("none", "PCA", "HillSmith", "dudi.mix"))
    stop("Argument ordination must be one of the following: 'none', 'PCA', 'HillSmith' or 'dudi.mix'")
  if (!is.null(noAxes) & ordination == "none")
    message(paste("You chose to use", noAxes, "ordination axes, but no ordination technique.\n",
                  "Ordination will be skipped"))
  if (HVmethod %in% c("box", "gaussian")) {
    if ((!is.null(bwHV1) | !is.null(bwHV2)) & freeBW) {
      message(paste("You provided bandwidth values, but freeBW is TRUE.\n",
                    "Bandwidths will be estimated per dimension using the default estimator.",
                    "See ?hypervolume::estimate_bandwidth"))
    }

    if (all(is.null(bwHV1), is.null(bwHV2))) {
      message("Only one set of bandwidths provided, both will be re-estimated.\n",
              "If this is not intended provide both sets of bandwidths.")
    }
  }

  if (is.null(init.vars)) {
    init.vars = c(1:ncol(HVdata1))
  } else{
    if (!all(is(init.vars, "numeric"), is(init.vars, "integer")))
      stop("init.vars must be a numeric/integer vector of columns indices")
  }

  if (!exists("HVidvar"))
    stop("Provide column number that contains HV IDs")

  if (!dir.exists(outputs.dir)) {
    dir.create(outputs.dir)
  }

  if (any(is(HVidvar, "integer"), is(HVidvar, "numeric"))) {
    HVnames <- unique(c(as.character(HVdata1[, HVidvar]), as.character(HVdata2[, HVidvar])))
    if (length(HVnames) > 2)
      stop("Only 2 hypervolumes can be used and HVidvar contains >2 IDs")
  } else
    stop("HVidvar should be the column number")

  if (length(HVidvar) > 1) stop("Provide one single column for HV IDs")

  ## check if variables contain the ID variable; if so remove it
  if (HVidvar %in% init.vars) init.vars = init.vars[-which(init.vars == HVidvar)]

  HVdata1$Type <- HVnames[1]
  HVdata2$Type <- HVnames[2]

  init.vars2 <- sort(c(init.vars, which(names(HVdata1) == "Type")))

  ## Joining tables and redo init.vars
  big.table <- rbind(HVdata1[, init.vars2], HVdata2[, init.vars2])
  init.vars <- grep("Type", names(big.table), invert = TRUE)

  if (any(is.na(big.table)))
    stop("NAs must be removed prior to computations")

  ## In case PFG relative abundances are to be used, there's no need to rescale.
  ## To avoid wrapping the whole function in the "if" the dataframe used in the PCAs is replaced
  if (do.scale) {
    big.table <- .scaleVars(big.table, init.vars)
  }

  ## ----------------------------------------------------------------------------
  ## ORDINATION only one ordination is calculated on the whole dataset
  if (ordination == "none") {
    HVpoints <- big.table[, init.vars]
    noAxes <- length(init.vars)

    if (noAxes == 1) stop("Only one dimension was selected. Please select at least two.")
    if (noAxes > 8) {
      warning("A minimum total explained variance of 90% is explained by >8 PCs.\n
                  Constraining no. of PCs to 8 anyway.")
      noAxes <- 8
    }
  } else {
    ordi.list <- HVordination(datatable = big.table, init.vars = init.vars, ordination = ordination,
                              HVidvar = which(names(big.table) == "Type"),
                              noAxes = noAxes, plotOrdi = plotOrdi, saveOrdi = saveOrdi,
                              outputs.dir = outputs.dir, file.suffix = paste(file.suffix, HVnames[1], HVnames[2], sep = "_"))

    HVpoints <- ordi.list[[1]]
    noAxes <- ordi.list[[2]]

    rm(ordi.list); gc(reset = TRUE)  ## clean workspace and memory
  }

  if (HVmethod %in% c("box", "gaussian")) {
    if (any(is.null(bwHV1), is.null(bwHV2)) | isTRUE(freeBW)) {
      ## if necessary estimate bandwidth for each HV (see Blonder et al. 2017)
      message("Bandwidth values will be calculated using the default estimator. See ?estimate_bandwidth")
      bwHV1 <- estimate_bandwidth(HVpoints[which(big.table$Type == HVnames[1]), 1:noAxes])
      bwHV2 <- estimate_bandwidth(HVpoints[which(big.table$Type == HVnames[2]), 1:noAxes])
    }
  }


  ## ----------------------------------------------------------------------------
  ## HYPERVOLUMES
  ## Do not parellelise - generates random NA's in the data for some reason
  for (i in 1:no.runs) {
    out <- .HVcalc(big.table, HVpoints, noAxes, HVmethod, ordination, HVnames, bwHV1, bwHV2, verbose, ...)
    HV1.disjfact <- out$HV1.disjfact
    HV2.disjfact <- out$HV2.disjfact

    if (freeBW == TRUE & HVmethod == "box") {
      while (HV1.disjfact >= 0.9 | HV2.disjfact >= 0.9) {
        bwHV1 <- bwHV1 + 0.05
        bwHV2 <- bwHV2 + 0.05
        out <- .HVcalc(big.table, HVpoints, noAxes, HVmethod, ordination, HVnames, bwHV1, bwHV2, verbose, ...)
        HV1.disjfact <- out$HV1.disjfact
        HV2.disjfact <- out$HV2.disjfact
      }
    }

    volumes <- data.frame("Dimensionality" = c(out$HV1@Dimensionality, out$HV2@Dimensionality),
                          "Volume" = c(out$HV1@Volume, out$HV2@Volume),
                          "Bandwidth" = out$Bandwidth,
                          "DisjunctFactor" = c(out$HV1.disjfact, out$HV2.disjfact),
                          "SVM_nu" = out$SVM_nu,
                          "SVM_gamma" = out$SVM_gamma,
                          row.names = c(HVnames[1], HVnames[2]))

    vol_comparison <- get_volume(out$volume.set)
    names(vol_comparison) = c(paste0("Volume_HV1_", HVnames[1]),
                              paste0("Volume_HV2_", HVnames[2]),
                              "Intersection",
                              "Union",
                              paste0("Unique_vol_HV1_", HVnames[1]),
                              paste0("Unique_vol_HV2_", HVnames[2]))

    vol_comparison <- data.frame(t(as.matrix(vol_comparison)),
                                 "MinDist" = out$hv.min.dist,
                                 "CentroidDist" = out$hv.centroid.dist)

    saveRDS(volumes, file = file.path(outputs.dir, paste0(file.suffix, "_HVdetails_", i,".rds")))
    saveRDS(vol_comparison, file = file.path(outputs.dir, paste0(file.suffix, "_Intersection_results_", i,".rds")))

    if (isTRUE(plotHV)) {
      cairo_pdf(filename = file.path(outputs.dir, paste0(file.suffix, "_Hypervolumes_", i,".pdf")),
                onefile = TRUE,
                width = 5, height = 5)
      plot(out$HV1)
      plot(out$HV2)
      plot(out$volume.set)
      graphics.off()
    }
  }
  paste("***done***")

}


#' CALCULATE ORDINATION FOR HYPERVOLUME COMPUTATION
#'
#' Calculates an ordination across two datasets prior to calculating two hypervolumes (one from each dataset)
#'
#' @param datatable is a \code{data.table} contained the raw data for the two hypervolumes
#'
#' @inheritParams hypervolumes
#'
#' @export
#'
#' @import data.table
#' @importFrom ade4 dudi.hillsmith dudi.mix
#' @importFrom stats prcomp predict
#'
#' @return a \code{list} of points used to build hypervolumes and final number of axes.

HVordination <- function(datatable, HVidvar, init.vars = NULL, ordination = "PCA",
                         noAxes = NULL, plotOrdi = TRUE, outputs.dir = NULL, file.suffix = NULL, saveOrdi = TRUE) {

  if (!ordination %in% c("PCA", "HillSmith", "dudi.mix")) stop("Argument ordination must be one of the following: 'none', 'PCA', 'HillSmith' or 'dudi.mix'")
  if (!is.null(noAxes) & ordination == "none")
    message(paste("You chose to use", noAxes, "ordination axes, but no ordination technique.\n
                  Ordination will be skipped"))
  if (saveOrdi) {
    if (is.null(outputs.dir) | is.null(file.suffix)) stop("Provide an output directory and/or file name")
  }

  if (!exists("HVidvar")) stop("Provide column number that contains HV IDs")
  if (class(HVidvar) == "integer" | class(HVidvar) == "numeric") {
    HVnames <- unique(datatable[, HVidvar])
    if (length(HVnames) > 2) stop("Only 2 hypervolumes can be used and HVidvar contains >2 IDs")
  } else stop("HVidvar should be the column number")

  if (length(HVidvar) > 1) stop("Provide one single column for HV IDs")


  if (is.null(init.vars)) {
    init.vars = c(1:ncol(datatable))
    message("Using all variables in data")
  } else{
    if (class(init.vars) != "numeric" & class(init.vars) != "integer") stop("init.vars must be a numeric/integer vector of columns indices")
  }

  ## check if variables contain the ID variable; if so remove it
  if (HVidvar %in% init.vars) init.vars = init.vars[-which(init.vars == HVidvar)]

  if (ordination == "PCA") {
    ## Q-mode; centered and using already scaled data
    ordi <- prcomp(datatable[, init.vars], center = TRUE, scale. = FALSE)
    fscores <- predict(ordi)

    ## correlations between variables and first axis
    # as.data.frame(ordi$rotation[order(abs(ordi$rotation[,1]), decreasing = TRUE),1])

    ## Extracting the % of explained variance per PC
    PEV <- ordi$sdev^2/sum(ordi$sdev^2)

    ## Extracting eigenvalues
    eigenv <- ordi$sdev^2
  }

  if (ordination == "HillSmith") {
    ordi <- dudi.hillsmith(datatable[, init.vars], scannf = FALSE, nf = length(init.vars))
    fscores <- ordi$li

    ## Extracting % of explained variable by PC
    PEV <- ordi$eig/sum(ordi$eig)

    ## Extracting eigenvalues
    eigenv <- ordi$eig
  }

  if (ordination == "dudi.mix") {
    ordi <- dudi.mix(datatable[, init.vars], scannf = FALSE, nf = length(init.vars))
    fscores <- ordi$li

    ## Extracting % of explained variable by PC
    PEV <- ordi$eig/sum(ordi$eig)

    ## Extracting eigenvalues
    eigenv <- ordi$eig
  }

  ## If the no of PCs is not determined, calculate how many are needed to have >=90% of the variance explained.
  if (is.null(noAxes)) {
    variance <- 0
    noAxes <- NA
    for (i in 1:length(PEV)) {
      variance <- variance + PEV[i]
      if (variance >= 0.9) {
        noAxes <- i
        break
      }
    }

    if (noAxes == 1) {
      warning("The first PC explains >= 90% of the total variance.\n
                  The two first PCs will be used as the dimensions for a 2D hypervolume nevertheless")
      noAxes = 2
    }

    if (noAxes > 8) {
      warning("A minimum total explained variance of 90% is explained by >8 PCs.\n
                  Constraining no. of PCs to 8.")
      noAxes <- 8
    }
  }

  ## PCA-related plots
  if (isTRUE(plotOrdi)) {
    .ordinationPlots(ordi, fscores, PEV, datatable, HVidvar, HVnames,
                     ordination, outputs.dir, file.suffix)
  }

  ## Saving PCA summary outputs
  if (saveOrdi) {
    if (ordination %in% c("HillSmith", "dudi.mix")) {

      summPCA <- rbind("Standard deviation" = sqrt(eigenv),
                       "Proportion of Variance" = PEV,
                       "Cumulative Proportion" = cumsum(PEV),
                       "Eigenvalues" = eigenv)
      saveRDS(summPCA, file.path(outputs.dir, paste(file.suffix, "OrdinationSumm.rds", sep = "_")))
    } else {
      summPCA <- rbind(summary(ordi)$importance, Eigenvalues = eigenv)
      saveRDS(summPCA, file.path(outputs.dir, paste(file.suffix, "OrdinationSumm.rds", sep = "_")))
    }
  }

  HVpoints <- fscores
  return(list("HVpoints" = HVpoints, "noAxes" = noAxes))
}

#' ORDINATION PLOT SAVING FUNCTION
#'
#' Calculates an ordination across two datasets prior to calculating two hypervolumes (one from each dataset)
#'
#' @param datatable is a \code{data.table} contained the raw data for the two hypervolumes
#' @param ordi the ordination object
#' @param fscores factor scores
#' @param PEV percent explained variance per PC
#' @param HVnames character vector of hypervolume names
#'
#' @inheritParams hypervolumes
#'
#' @import data.table
#' @importFrom ade4 scatter
#' @importFrom stats biplot
#' @importFrom graphics barplot layout par points
#' @importFrom grDevices cairo_pdf graphics.off
#'
#' @return NULL, just saves plots to a PDF

.ordinationPlots <- function(ordi, fscores, PEV, datatable, HVidvar, HVnames,
                             ordination, outputs.dir, file.suffix) {
  plotFUN <- if (ordination %in% c("PCA", "HillSmith")) {
    biplot
  } else {
    scatter
  }

  plotTitle <- paste(HVnames[1], "(black) and", HVnames[2], "(red)")
  colVect <- ifelse(datatable[, HVidvar] == HVnames[1], "black", "red")
  pchVect <- ifelse(datatable[, HVidvar] == HVnames[1], 19, 1)

  cairo_pdf(filename = file.path(outputs.dir, paste(file.suffix, "Ordination.pdf", sep = "_")), width = 10, height = 10)
  layout(matrix(c(1:4), nrow = 2, ncol = 2, byrow = TRUE))
  sets <- options(warn = -1)  ## suppressing warnings about zero-length arrows
  plotFUN(ordi, choices = c(1,2))
  plotFUN(ordi, choices = c(2,3))
  plot(fscores[, 1], fscores[, 2],
       xlim = range(fscores[,1]), ylim = range(fscores[,2]),
       col = colVect, pch = pchVect,
       main = plotTitle, xlab = "PC1", ylab = "PC2")
  options(sets)
  par(new = FALSE)
  barplot(PEV,
          main = paste(file.suffix, "OrdiScreeplot.pdf", sep = "_"),
          ylab = "% explained variance", xlab = "Principal components")
  graphics.off()
}

#' HYPERVOLUME COMPUTATION WRAPPER
#'
#' Calculates two hypervolumes and compares them in terms of
#'   overlap and distance. When appropriate (i.e. \code{HVmethod == "box"}),
#'   it also calculates the disjunct factor and repeats the calculations
#'   with larger bandwidths (if \code{freeBW == TRUE})
#'
#' @param big.table is a \code{data.table} contained the raw data for the two hypervolumes.
#'   Used to subset \code{HVpoints} according to HVnames, so rows have to correspond
#' @param HVpoints \code{data.table} of points used to build hypervolumes and final number of axes.
#' @param HVnames \code{character} vector of length 2 of hypervolume names.
#' @inheritParams hypervolumes
#'
#' @import data.table
#' @importFrom hypervolume hypervolume hypervolume_distance get_volume estimate_bandwidth hypervolume_set
#'
#' @return a \code{list} with entries:
#'  \itemize{
#'    \item 'HV1'
#'    \item 'HV2'
#'    \item 'volume.set' - the result of \code{hypervolume::hypervolume_set}
#'    \item 'hv.centroid.dist' and 'hv.min.dist' - the result of \code{hypervolume::hypervolume_distance(..., type = "centroid")} and \code{hypervolume::hypervolume_distance(..., type = "minimum")}
#'    \item 'Bandwidth' - \code{NA} if \code{HVmethod} is not "box" or "gaussian")
#'    \item 'HV1.disjfact' and 'HV2.disjfact' (\code{NA} if \code{HVmethod} is not "box")
#'    \item 'SVM_nu' and 'SVM_gama' (\code{NA} if \code{HVmethod} is "box" or "gaussian")
#'  }

.HVcalc <- function(big.table, HVpoints, noAxes, HVmethod, ordination,
                    HVnames, bwHV1, bwHV2, verbose, ...) {
  if (ordination != "none") {
    HV1name <- paste(HVnames[1], "- Ordination factor scores")
    HV2name <- paste(HVnames[2], "- Ordination factor scores")
  } else {
    HV1name <- HVnames[1]
    HV2name <- HVnames[2]
  }

  if (HVmethod %in% c("box", "gaussian")) {
    HV1name <- paste(HV1name,"- fixed bandwidth to", bwHV1)
    HV2name <- paste(HV2name,"- fixed bandwidth to", bwHV2)

    HV1 <- hypervolume(data = HVpoints[which(big.table$Type == HVnames[1]), 1:noAxes],
                       method = HVmethod, kde.bandwidth = bwHV1,
                       name = HV1name,
                       verbose = verbose, ...)

    HV2 <- hypervolume(data = HVpoints[which(big.table$Type == HVnames[2]), 1:noAxes],
                       method = HVmethod, kde.bandwidth = bwHV2,
                       name = HV2name,
                       verbose = verbose, ...)
  } else {
    HV1 <- hypervolume(data = HVpoints[which(big.table$Type == HVnames[1]), 1:noAxes],
                       method = HVmethod, name = HV1name,
                       verbose = verbose, ...)

    HV2 <- hypervolume(data = HVpoints[which(big.table$Type == HVnames[2]), 1:noAxes],
                       method = HVmethod, name = HV2name,
                       verbose = verbose, ...)
  }

  volume.set <- hypervolume_set(HV1, HV2, check.memory = FALSE, verbose = verbose)

  hv.centroid.dist <- hypervolume_distance(HV1, HV2, type = "centroid")
  hv.min.dist <- hypervolume_distance(HV1, HV2, type = "minimum", check.memory = FALSE)

  if (HVmethod %in% c("box", "gaussian")) {
    Bandwidth <- c(HV1@Parameters$kde.bandwidth[1], HV2@Parameters$kde.bandwidth[1])
    SVM_nu <- NA
    SVM_gamma <- NA
  } else {
    Bandwidth <- NA
    SVM_nu <- HV1@Parameters$svm.nu
    SVM_gamma <- HV1@Parameters$svm.gamma
  }

  ## Disjunct factor is no longer an output (hypervolume v 2.0.8) -
  ## now calculated "manually" when using the box method (not relevant in gaussian/svm)
  if (HVmethod == "box") {
    HV1.disjfact <- HV1@Volume / (nrow(HVpoints[which(big.table$Type == HVnames[1]), 1:noAxes]) * prod(2 * HV1@Parameters$kde.bandwidth))
    HV2.disjfact <- HV2@Volume / (nrow(HVpoints[which(big.table$Type == HVnames[2]), 1:noAxes]) * prod(2 * HV2@Parameters$kde.bandwidth))
  } else {
    HV1.disjfact <- NA
    HV2.disjfact <- NA
  }

  list("HV1" = HV1,
       "HV2" = HV2,
       "volume.set" = volume.set,
       "hv.centroid.dist" = hv.centroid.dist,
       "hv.min.dist" = hv.min.dist,
       "Bandwidth" = Bandwidth,
       "HV1.disjfact" = HV1.disjfact,
       "HV2.disjfact" = HV2.disjfact,
       "SVM_nu" = SVM_nu,
       "SVM_gamma" = SVM_gamma)
}


#' SCALE VARIABLES
#'
#' Scales variables/columns in \code{init.vars} using \code{scale(..., center = FALSE, scale = TRUE)}
#' Variables that only have 0s, will be assigned 1e-6
#'
#' @param datatable a table object coercible to \code{data.frame}
#' @inheritParams HVordination
#'
#' @return scaled data.frame
#'
.scaleVars <- function(datatable, init.vars) {
  datatable <- as.data.frame(datatable)
  ## Attributing a very small value (1e-6) to variables that have only 0s enables scaling of all variables w/o producing NaNs
  datatable[, which(colSums(datatable[, init.vars]) == 0)] <- 1e-6

  ## Scaling using root mean squares (by having center=FALSE, scale=TRUE) to avoid producing NaNs in constant variables
  ## Scaling needs to be done with both datasets together, otherwise intersections are forced
  datatable[, init.vars] <- as.data.frame(scale(datatable[, init.vars],
                                                center = FALSE, scale = TRUE))
  return(datatable)
}
