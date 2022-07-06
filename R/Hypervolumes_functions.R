globalVariables(c("HVidvar", "..init.vars", "..init.vars2"))

#' CALCULATE AND COMPARE TWO HYPERVOLUMES
#'
#' Calculates two \code{hypervolumes} on two sets of "raw" data or on the factor scores of
#' an ordination done on the entire dataset (across the two sets of data).
#'
#' @param HVdata1 the \code{data.table} (or a object that can be converted to \code{data.table}) from which the first hypervolume will be calculated
#' @param HVdata2 the \code{data.table} (or a object that can be converted to \code{data.table}) from  which the second hypervolume will be calculated.
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
#' @param saveOrdi activates/deactivates saving ordination object
#' @param saveOrdiSumm activates/deactivates saving ordination summary outputs
#' @param outputs.dir is the directory to store results
#' @param file.suffix is a character string used as a suffix in file names
#' @param verbose is passed to \code{hypervolume::*} functions to control printing of diagnostic outputs
#' @param plotHVDots a \code{list} of further arguments passed to \code{hypervolume:::plot.HypervolumeList}.
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
#' @importFrom hypervolume plot.HypervolumeList
#' @importFrom methods is

hypervolumes <- function(HVdata1, HVdata2, HVidvar, ordination = "PCA", init.vars = NULL,
                         noAxes = NULL, do.scale = FALSE, HVmethod = "box",
                         freeBW = FALSE, bwHV1 = NULL, bwHV2 = NULL,
                         no.runs = 1, plotOrdi = TRUE, plotHV = TRUE, saveOrdi = FALSE, saveOrdiSumm = TRUE,
                         outputs.dir, file.suffix, verbose = TRUE, plotHVDots = NULL, ...) {
  ## do some checks:
  if (!is(HVdata1, "data.table")) {
    HVdata1 <- as.data.table(HVdata1)
  }
  if (!is(HVdata2, "data.table")) {
    HVdata2 <- as.data.table(HVdata2)
  }
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

    if (any(is.null(bwHV1), is.null(bwHV2))) {
      message("One or no sets of bandwidths provided; both will be re-estimated.\n",
              "If this is not intended provide both sets of bandwidths.")
    }
  }

  if (is.null(init.vars)) {
    init.vars <- c(1:ncol(HVdata1))
  } else{
    if (!all(is(init.vars, "numeric"), is(init.vars, "integer")))
      stop("init.vars must be a numeric/integer vector of columns indices")
  }

  if (!exists("HVidvar"))
    stop("Provide column number that contains HV IDs")

  if (!dir.exists(outputs.dir)) {
    dir.create(outputs.dir)
  }

  if (length(HVidvar) > 1) stop("Provide one single column for HV IDs")

  if (any(is(HVidvar, "integer"), is(HVidvar, "numeric"))) {
    HVnames <- unique(c(as.character(HVdata1[[HVidvar]]), as.character(HVdata2[[HVidvar]])))
    if (length(HVnames) > 2)
      stop("Only 2 hypervolumes can be used and HVidvar contains >2 IDs")
  } else
    stop("HVidvar should be the column number")

  ## check if variables contain the ID variable; if so remove it
  if (HVidvar %in% init.vars) init.vars <- init.vars[-which(init.vars == HVidvar)]

  HVdata1$Type <- HVnames[1]
  HVdata2$Type <- HVnames[2]

  init.vars2 <- sort(c(init.vars, which(names(HVdata1) == "Type")))

  ## Joining tables and redo init.vars
  big.table <- rbind(HVdata1[, ..init.vars2], HVdata2[, ..init.vars2])
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
    HVpoints <- big.table[, ..init.vars]
    noAxes <- length(init.vars)

    if (noAxes == 1) stop("Only one dimension was selected. Please select at least two.")
    if (noAxes > 8) {
      warning("Hypervolumes with > 8 dimensions are not allowed\n
                  Constraining no. of PCs to 8")
      noAxes <- 8
    }
  } else {
    ordi.list <- HVordination(datatable = big.table, init.vars = init.vars, ordination = ordination,
                              HVidvar = which(names(big.table) == "Type"), noAxes = noAxes,
                              plotOrdi = plotOrdi, saveOrdi = saveOrdi, saveOrdiSumm = saveOrdiSumm,
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
      do.call(plot, args = c(list(x = out$HV1), plotHVDots))
      do.call(plot, args = c(list(x = out$HV2), plotHVDots))
      do.call(plot, args = c(list(x = out$volume.set), plotHVDots))
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
#' @return a \code{list} of points used to build hypervolumes and final number of axes.
#'
#' @export
#'
#' @import data.table
#' @importFrom ade4 dudi.hillsmith dudi.mix
#' @importFrom stats prcomp predict
#' @importFrom methods is

HVordination <- function(datatable, HVidvar, init.vars = NULL, ordination = "PCA",
                         noAxes = NULL, plotOrdi = TRUE, outputs.dir = NULL,
                         file.suffix = NULL, saveOrdi = FALSE, saveOrdiSumm = TRUE) {

  if (!is(datatable, "data.table")) {
    datatable <- as.data.table(datatable)
  }

  if (!ordination %in% c("PCA", "HillSmith", "dudi.mix")) stop("Argument ordination must be one of the following: 'none', 'PCA', 'HillSmith' or 'dudi.mix'")
  if (!is.null(noAxes) & ordination == "none")
    message(paste("You chose to use", noAxes, "ordination axes, but no ordination technique.\n
                  Ordination will be skipped"))
  if (saveOrdi | saveOrdiSumm) {
    if (is.null(outputs.dir) | is.null(file.suffix)) stop("Provide an output directory and/or file name")
  }

  if (length(HVidvar) > 1) stop("Provide one single column for HV IDs")

  if (!exists("HVidvar")) stop("Provide column number that contains HV IDs")

  if (any(is(HVidvar, "integer"), is(HVidvar, "numeric"))) {
    HVnames <- unique(datatable[[HVidvar]])
    if (length(HVnames) > 2) warning("Only 2 hypervolumes can be compared and HVidvar contains >2 IDs!")
  } else stop("HVidvar should be the column number")

  if (is.null(init.vars)) {
    init.vars = c(1:ncol(datatable))
    message("Using all variables in data")
  } else{
    if (class(init.vars) != "numeric" & class(init.vars) != "integer") stop("init.vars must be a numeric/integer vector of columns indices")
  }

  ## check if variables contain the ID variable; if so remove it
  if (HVidvar %in% init.vars) init.vars <- init.vars[-which(init.vars == HVidvar)]

  if (ordination == "PCA") {
    ## Q-mode; centered and using already scaled data
    ordi <- prcomp(datatable[, ..init.vars], center = TRUE, scale. = FALSE)
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

  ## Saving PCA object
  if (saveOrdi) {
    saveRDS(ordi, file.path(outputs.dir, paste(file.suffix, "OrdinationObj.rds", sep = "_")))
  }

  ## Saving PCA summary outputs
  if (saveOrdiSumm) {
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
  return(as.data.table(datatable))
}


#' PLOT 3D HYPERVOLUMES WITH ORDINATION LOADINGS AND POST-HOC FITTED VECTORS
#'
#' Plots two 3D hypervolumes that are based on ordination scores, showing
#' variable loadings and vectors fitted post-hoc to the ordination
#'
#' @param HVlist a \code{HypervolumeList}
#' @param loadings_coords a data.table or data.frame of ordination factor
#'   loadings for the 3 axes to be plotted
#' @param PHvect_coords a data.table or data.frame of post-hoc fitted vectors
#'   loadings for the 3 axes to be plotted
#' @param loadings_labels character vector of labels for each ordination factor.
#'   Defaults to \code{row.names(loadings_coords)}
#' @param PHvect_labels character vector of labels for each post-hoc fitted vector.
#'   Defaults to \code{row.names(PHvect_coords)}
#' @param loadings_col a colour for ordination factor arrows and labels.
#' @param PHvect_col a colour for post-hoc fitted vector arrows and labels.
#' @param cex.label passed to \code{basicPlotteR::addTextLabels}
#' @param lwd passed to \code{scatterplot3d::scatterplot3d}
#' @param ... arguments passed to \code{.plot.HypervolumeListEdited}
#'
#' @importFrom basicPlotteR addTextLabels
#' @export

plotHypervolumes3D <- function(HVlist, loadings_coords = NULL, PHvect_coords = NULL,
                               loadings_labels = NULL, PHvect_labels = NULL,
                               loadings_col = "grey20", cex.label = 1, lwd = 1,
                               PHvect_col = "darkgreen", ...) {

  if (!is.null(loadings_coords)) {
    if (ncol(loadings_coords) != 3) {
      stop("loadings_coords must have three columns for 3D plotting")
    }
    if (is.null(loadings_labels)) {
      loadings_labels <- paste0("loadingVar_", row.names(loadings_coords))
    }
  }

  if (!is.null(PHvect_coords)){
    if (ncol(PHvect_coords) != 3) {
      stop("PHvect_coords must have three columns for 3D plotting")
    }
    if (is.null(PHvect_labels)) {
      PHvect_labels <-  paste0("PHVec_", row.names(PHvect_coords))
    }
  }

  HVPlot <- .plot.HypervolumeListEdited(x = HVlist, show.3d = TRUE, rglGraphics = FALSE, ...)

  if (!is.null(loadings_coords)) {
    loadings_coords$label <- loadings_labels
    loadings_coords$type <- "loadings"
    ## adding factor loadings
    HVPlot$points3d(loadings_coords[, 1:3],
                    type = "l",
                    col = makeTransparent(loadings_col, alpha = 200), lwd = lwd)
  } else {
    loadings_coords <- data.table(matrix(data = "", ncol = HVlist@HVList[[1]]@Dimensionality + 2,
                                         dimnames = list(rownames = NULL,
                                                         colnames = c(colnames(HVlist@Data), "label", "type"))))[0]
  }

  if (!is.null(PHvect_coords)) {
    ## create a label column
    PHvect_coords$label <- PHvect_labels
    PHvect_coords$type <- "PHvect"

    ## adding trait vectors
    HVPlot$points3d(PHvect_coords[, 1:3],
                    type = "l",
                    col = makeTransparent(PHvect_col, alpha = 200), lwd = lwd)
  } else {
    PHvect_coords <- data.table(matrix(data = "", ncol = HVlist@HVList[[1]]@Dimensionality + 2,
                                       dimnames = list(rownames = NULL,
                                                       colnames = c(colnames(HVlist@Data), "label", "type"))))[0]
  }

  ## plotting labels
  ## labels need to be plotted altogether to avoid overlapping
  ## function doesn't take a colour vector, so the labels that have another colour are left empty
  allVectors <- rbind(loadings_coords, PHvect_coords)

  if (nrow(allVectors)) {
    ## text coordinates
    textCoords <- HVPlot$xyz.convert(allVectors[rowSums(allVectors[,1:3]) != 0, 1:3])

    ## factor loading labels
    plotlabels <- allVectors$label[rowSums(allVectors[,1:3]) != 0]
    plotlabels[allVectors$type[rowSums(allVectors[,1:3]) != 0] == "PHvect"] <- ""
    addTextLabels(xCoords = textCoords$x, yCoords = textCoords$y,
                  labels = plotlabels, cex.label = cex.label, col.label = loadings_col, col.line = loadings_col, lwd = 0.15, lty = 2,
                  avoidPoints = TRUE)

    ## trait labels
    plotlabels <- allVectors$label[rowSums(allVectors[,1:3]) != 0]
    plotlabels[allVectors$type[rowSums(allVectors[,1:3]) != 0] == "loadings"] <- ""
    addTextLabels(xCoords = textCoords$x, yCoords = textCoords$y,
                  labels = plotlabels, cex.label = cex.label,
                  col.label = PHvect_col, col.line = PHvect_col, lwd = 0.15, lty = 2,
                  avoidPoints = TRUE)
  }
}

#' Edited version of \code{hypervolume::plot.HypervolumeList}
#'
#' Added centroid colours, corrected cex.data
#' graphics changed to \code{scatterplot3d}, instead of \code{rgl}

#' @inheritParams hypervolume::plot.HypervolumeList
#' @param rglGraphics if TRUE will plot 3D hypervolumes using \code{rgl}.
#'   Otherwise, \code{scatterplot3d} is used
#' @param centroid.cols a vector of as many colours as hypervolumes used
#'   to colour centroids
#' @param ... passed to \code{points} if \code{show.3d == FALSE}, or
#'   otherwise to \code{scatterplot3d}.
#'
#' @importFrom grDevices col2rgb rainbow rgb
#' @importFrom graphics axis box contour legend mtext text
#' @importFrom stats quantile
#' @importFrom scatterplot3d scatterplot3d
#' @importFrom rgl plot3d points3d mtext3d
#' @importFrom MASS kde2d
#' @export

.plot.HypervolumeListEdited <- function(x, show.3d = FALSE, rglGraphics = TRUE, plot.3d.axes.id = NULL, show.axes = TRUE,
                                        show.frame = TRUE, show.random = TRUE, show.density = TRUE,
                                        show.data = TRUE, names = NULL, show.legend = TRUE, limits = NULL,
                                        show.contour = TRUE, contour.lwd = 1.5, contour.type = "kde",
                                        contour.alphahull.alpha = 0.25, contour.ball.radius.factor = 1,
                                        contour.kde.level = 0.01, contour.raster.resolution = 100,
                                        show.centroid = TRUE, cex.centroid = 2,
                                        colors = rainbow(floor(length(x@HVList) * 1.5), alpha = 0.8),
                                        centroid.cols = rainbow(floor(length(x@HVList) * 1.5), alpha = 0.8),
                                        point.alpha.min = 0.2, point.dark.factor = 0.5,
                                        cex.random = 0.5, cex.data = 0.75, cex.axis = 0.75, cex.names = 1,
                                        cex.legend = 0.75, num.points.max.data = 1000, num.points.max.random = 2000,
                                        reshuffle = TRUE, plot.function.additional = NULL, verbose = FALSE,
                                        ...) {
  sapply(x@HVList, function(z) {
    if (verbose == TRUE) {
      cat(sprintf("Showing %d random points of %d for %s\n",
                  min(nrow(z@RandomPoints), num.points.max.random),
                  nrow(z@RandomPoints), z@Name))
    }
    if (show.data && length(z@Data) > 0) {
      npd <- ifelse(all(is.nan(z@Data)), 0, nrow(z@Data))
      if (verbose == TRUE) {
        cat(sprintf("Showing %d data points of %d for %s\n",
                    min(num.points.max.data, npd), npd, z@Name))
      }
    }
  })
  if (!requireNamespace("alphahull", quietly = TRUE)) {
    warning("The package 'alphahull' is needed for contour plotting with contour.type='alphahull'. Please install it to continue.\n\n *** Temporarily setting contour.type='kde'.",
            call. = FALSE)
    contour.type <- "kde"
  }
  alldims = sapply(x@HVList, function(z) {
    z@Dimensionality
  })
  allnames = sapply(x@HVList, function(z) {
    z@Name
  })
  stopifnot(all(alldims[1] == alldims))
  all <- NULL
  alldata <- NULL
  for (i in 1:length(x@HVList)) {
    ivals = sample(nrow(x@HVList[[i]]@RandomPoints), min(c(num.points.max.random,
                                                           nrow(x@HVList[[i]]@RandomPoints))))
    subsampledpoints = data.frame(x@HVList[[i]]@RandomPoints[ivals,
                                                             , drop = FALSE])
    densityvals = x@HVList[[i]]@ValueAtRandomPoints[ivals]
    if (nrow(subsampledpoints) > 0) {
      subsampledpoints = cbind(subsampledpoints, ID = rep(i,
                                                          nrow(subsampledpoints)), Density = (densityvals -
                                                                                                min(densityvals, na.rm = TRUE))/(max(densityvals,
                                                                                                                                     na.rm = TRUE) - min(densityvals, na.rm = TRUE)))
      subsampledpoints[is.nan(subsampledpoints[, "Density"]),
                       "Density"] <- 1
      all <- rbind(all, subsampledpoints)
    }
    thisdata = x@HVList[[i]]@Data
    alldata <- rbind(alldata, cbind(thisdata, ID = rep(i,
                                                       nrow(thisdata))))
  }
  alldata <- as.data.frame(alldata)
  if (num.points.max.data < nrow(alldata) && !is.null(num.points.max.data)) {
    alldata <- alldata[sample(nrow(alldata), min(c(num.points.max.data,
                                                   nrow(alldata)))), ]
  }
  if (is.null(all)) {
    warning("No random points to plot.")
    if (is.null(dimnames(x@HVList[[1]]@RandomPoints)[[2]])) {
      all <- matrix(0, ncol = 2 + alldims, nrow = 1, dimnames = list(NULL,
                                                                     c(paste("X", 1:alldims, sep = ""),
                                                                       "ID", "Density")))
    }
    else {
      all <- matrix(0, ncol = 2 + alldims, nrow = 1, dimnames = list(NULL,
                                                                     c(dimnames(x@HVList[[1]]@RandomPoints)[[2]],
                                                                       "ID", "Density")))
    }
    all <- as.data.frame(all)
  }
  if (reshuffle == TRUE) {
    all <- all[sample(nrow(all), replace = FALSE), , drop = FALSE]
    alldata <- alldata[sample(nrow(alldata), replace = FALSE),
                       , drop = FALSE]
  }
  no_names_supplied = FALSE
  if (is.null(names)) {
    dn = dimnames(all)[[2]]
    names = dn[1:(ncol(all) - 2)]
    no_names_supplied = TRUE
  }
  if (!is.null(limits) & !is.list(limits)) {
    varlimlist = vector("list", ncol(all) - 2)
    for (i in 1:length(varlimlist)) {
      varlimlist[[i]] <- limits
    }
    limits = varlimlist
  }
  colorlist <- colors[all$ID]
  alphavals <- (all$Density - quantile(all$Density, 0.025,
                                       na.rm = T))/(quantile(all$Density, 0.975, na.rm = T) -
                                                      quantile(all$Density, 0.025, na.rm = T))
  alphavals[is.nan(alphavals)] <- 0.5
  alphavals[alphavals < 0] <- 0
  alphavals[alphavals > 1] <- 1
  alphavals <- point.alpha.min + (1 - point.alpha.min) * alphavals
  if (show.density == FALSE) {
    alphavals <- rep(1, length(colorlist))
  }
  for (i in 1:length(colorlist)) {
    colorlist[i] <- hypervolume:::rgb_2_rgba(colorlist[i], alphavals[i])
  }
  colorlistdata = colors[alldata$ID]
  for (i in 1:length(colorlistdata)) {
    colorlistdata[i] <- hypervolume:::rgb_2_set_hsv(colorlistdata[i], v = 1 -
                                                      point.dark.factor)
  }
  if (ncol(all) < 2) {
    stop("Plotting only available in n>=2 dimensions.")
  }
  if (show.3d == FALSE) {
    op = par(no.readonly = T)
    par(mfrow = c(ncol(all) - 2, ncol(all) - 2))
    par(mar = c(0, 0, 0, 0))
    par(oma = c(0.5, 0.5, 0.5, 0.5))
    for (i in 1:(ncol(all) - 2)) {
      for (j in 1:(ncol(all) - 2)) {
        if (j > i) {
          plot(all[, j], all[, i], type = "n",
               axes = FALSE, xlim = limits[[j]], ylim = limits[[i]],
               bty = "n")
          if (show.random == TRUE) {
            points(all[, j], all[, i], col = colorlist,
                   cex = cex.random, pch = list(...)$pch)
          }
          if (show.data & nrow(alldata) > 0) {
            points(alldata[, j], alldata[, i], col = colorlistdata,
                   cex = cex.data, pch = list(...)$pch)
          }
          if (show.centroid == TRUE) {
            for (whichid in 1:length(unique(all$ID))) {
              allss <- subset(all, all$ID == whichid)
              centroid_x <- mean(allss[, j], na.rm = TRUE)
              centroid_y <- mean(allss[, i], na.rm = TRUE)
              points(centroid_x, centroid_y, col = colors[whichid],
                     cex = cex.centroid, pch = list(...)$pch)
              points(centroid_x, centroid_y, col = "white",
                     cex = cex.centroid, pch = 1, lwd = 1.5)
            }
          }
          if (show.contour == TRUE) {
            for (whichid in 1:length(unique(all$ID))) {
              allss <- subset(all, all$ID == whichid)
              if (nrow(allss) > 0) {
                contourx <- allss[, j]
                contoury <- allss[, i]
                rp = cbind(contourx, contoury)
                vol_this = x@HVList[[whichid]]@Volume
                density_this = nrow(rp)/vol_this
                dim_this = x@HVList[[whichid]]@Dimensionality
                radius_critical <- density_this^(-1/dim_this) *
                  contour.ball.radius.factor
                if (contour.type == "alphahull") {
                  poly_outline = hypervolume:::do_outline_alpha(rp = rp,
                                                                alpha = contour.alphahull.alpha)
                  plot(poly_outline, add = TRUE, wpoints = FALSE,
                       wlines = "none", lwd = contour.lwd,
                       col = colors[whichid])
                }
                else if (contour.type == "ball") {
                  poly_outline <- hypervolume:::do_outline_ball(rp = rp,
                                                                radius = radius_critical)
                  sp::plot(poly_outline, add = TRUE,
                           lwd = contour.lwd, col = colors[whichid])
                }
                else if (contour.type == "kde") {
                  m_kde = kde2d(rp[, 1], rp[, 2], n = 50,
                                h = radius_critical)
                  contour(m_kde, add = TRUE, levels = contour.kde.level,
                          drawlabels = FALSE, lwd = contour.lwd,
                          col = colors[whichid])
                }
                else if (contour.type == "raster") {
                  poly_raster <- hypervolume:::do_outline_raster(as.matrix(rp),
                                                                 res = contour.raster.resolution)
                  sp::plot(poly_raster, add = TRUE, lwd = contour.lwd,
                           col = colors[whichid])
                }
              }
            }
          }
          if (!is.null(plot.function.additional)) {
            plot.function.additional(j, i)
          }
          if (show.frame == TRUE) {
            box()
          }
        }
        else if (j == i) {
          plot(0, 0, type = "n", xlim = c(0, 1),
               ylim = c(0, 1), axes = FALSE)
          text(0.5, 0.5, names[j], cex = cex.names)
        }
        else if (j == 1 & i == (ncol(all) - 2)) {
          plot(0, 0, type = "n", xlim = c(0, 1),
               ylim = c(0, 1), axes = FALSE)
          if (show.legend == TRUE) {
            legend("topleft", legend = allnames,
                   text.col = colors, bty = "n", cex = cex.legend)
          }
        }
        else {
          plot(0, 0, type = "n", axes = FALSE)
        }
        if (j == i + 1) {
          if (show.axes == TRUE) {
            axis(side = 1, cex.axis = cex.axis)
            axis(side = 2, cex.axis = cex.axis)
          }
        }
      }
    }
    par(op)
  }
  else {
    if (is.null(plot.3d.axes.id)) {
      plot.3d.axes.id = 1:3
    }
    if (no_names_supplied == TRUE) {
      axesnames <- names[plot.3d.axes.id]
    }
    else {
      axesnames <- names
    }
    if (length(plot.3d.axes.id) != 3) {
      stop("Must specify three axes")
    }
    if (show.density == TRUE) {
      for (i in 1:length(colorlist)) {
        colorlist[i] <- hypervolume:::rgb_2_set_hsv(colorlist[i], s = (alphavals[i]^2))
      }
    }
    if (rglGraphics) {
      plot3d(all[, plot.3d.axes.id], col = colorlist,
             expand = 1.05, xlab = axesnames[1], ylab = axesnames[2],
             zlab = axesnames[3], xlim = limits[[1]], ylim = limits[[2]],
             zlim = limits[[3]], size = cex.random, type = "p",
             box = show.frame, axes = show.axes)
    } else {
      s3d <- scatterplot3d(x = all[, plot.3d.axes.id][, 1], y = all[, plot.3d.axes.id][, 2], z = all[, plot.3d.axes.id][, 3],
                           color = colorlist, xlab=axesnames[1], ylab=axesnames[2], zlab=axesnames[3],
                           xlim=limits[[1]], ylim=limits[[2]], zlim=limits[[3]],
                           cex.symbols = cex.random, cex.axis = cex.axis, type = "p", ...)
    }

    if (show.legend == TRUE) {
      if (rglGraphics) {
        for (i in 1:length(allnames)) {
          mtext3d(allnames[i], edge = "x-+",
                  line = 1 + i * cex.legend * 1.25, color = colors[i],
                  cex = cex.legend)
        }
      } else {
        for (i in 1:length(allnames)) {
          mtext(allnames[i], side = 1,line=1+i*cex.legend*1.25, col=colors[i], cex=cex.legend)
        }
      }
    }

    if (show.data) {
      if (!any(is.nan(as.matrix(alldata[, plot.3d.axes.id])))) {
        if (rglGraphics) {
          points3d(x = alldata[, plot.3d.axes.id[1]],
                   y = alldata[, plot.3d.axes.id[2]], z = alldata[, plot.3d.axes.id[3]],
                   col = colorlistdata,
                   cex = cex.data, pch = list(...)$pch)
        } else {
          s3d$points3d(x = alldata[,plot.3d.axes.id[1]], y = alldata[,plot.3d.axes.id[2]], z = alldata[,plot.3d.axes.id[3]],
                       col = colorlistdata, cex = cex.data, pch = list(...)$pch)
        }
      }

    }
    if (show.centroid == TRUE) {
      for (whichid in 1:length(unique(all$ID))) {
        allss <- subset(all, all$ID == whichid)
        centroid_1 <- mean(allss[, plot.3d.axes.id[1]],
                           na.rm = TRUE)
        centroid_2 <- mean(allss[, plot.3d.axes.id[2]],
                           na.rm = TRUE)
        centroid_3 <- mean(allss[, plot.3d.axes.id[3]],
                           na.rm = TRUE)
        if (rglGraphics) {
          points3d(x = centroid_1, y = centroid_2,
                   z = centroid_3, col = centroid.cols[whichid], size = cex.centroid)
        } else {
          s3d$points3d(x = centroid_1, y = centroid_2, z = centroid_3,
                       col = centroid.cols[whichid], cex = cex.centroid, pch = list(...)$pch)
          s3d$points3d(x = centroid_1, y = centroid_2, z = centroid_3,
                       col = "white", cex = cex.centroid, pch = 1, lwd = 1.5)
        }

      }
    }
  }
  if (!rglGraphics) {
    return(s3d)
  }
}
