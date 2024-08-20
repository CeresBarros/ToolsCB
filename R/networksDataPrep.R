##  ----------------------------------------------
##  NETWORKS DATA PROCESSING FUNCTIONS
##  ----------------------------------------------

#' IDENTIFY SPP BY CODE
#'
#' Function to Identify a spp by the code.
#' Author: João Braga 2016
#'
#' @template SPPCODE
#' @param SPPNAME Latin species name
#' @importFrom utils read.table
#' @export
whois <- function(SPPCODE = NULL, SPPNAME = NULL) {

  if (is.null(SPPCODE) & is.null(SPPNAME)) {
    stop("Must specify a species code or name(Genus_species)")
  }
  if (!is.null(SPPCODE) & !is.null(SPPNAME)) {
    stop("Must specify a species code or name(Genus_species)")
  }

  SppID <- read.table(file = "SppID.txt", header = TRUE)

  if (length(SPPCODE) > 1) {
    SPPCODE <- paste0(SPPCODE, "$", collapse = "|")
  }
  if (length(SPPNAME) > 1) {
    SPPNAME <- paste0(SPPNAME, "$", collapse = "|")
  }

  if (!is.null(SPPCODE)) {
    who <- SppID[grepl(pattern = SPPCODE,x = SppID$ID),]
  }
  if (!is.null(SPPNAME)){
    who <- SppID[grepl(pattern = SPPNAME,x = SppID$SPPname),]
  }

  return(who)
}

#' READ SEVERAL RDATA FILES TO A LIST
#'
#' Author: João Braga 2016
#'
#' @param file.list list of RDATA files
#' @export
fun.load.many.rdata <- function(file.list = NULL) {
  # function to read and assign rdata to a list
  if (is.null(file.list)) {
    stop("Must provide list of full path files")
  }
  e1 <- new.env()
  invisible(lapply(file.list, load, envir = e1))
  my_list = as.list(e1)

  return(my_list)
}


#' @title DETECT SPECIES PRESENCES IN PIXEL BY A THRESHOLD NUMBER OF PRESENCES
#'
#' @param species.dist a species distribution dbf file
#' @param threshold cell habitat percentage threshold to consider a species present
#' @param opt.only if `TRUE` only presences in optimal habitats are
#'   considered, if `FALSE` both secondary and optimal habitat presences are used
#'
#' @export
fun.PRESENCE.ABSENCE <- function(species.dist = NULL, threshold = NULL, opt.only = NULL) {
  if (is.null(species.dist)) stop("Must specify a species distribution dbf file")
  if (is.null(threshold)) stop("Must specify a cell habitat % threshold to consider a species present")
  if (is.null(opt.only)) stop("Must specify whether to use only optimal habitats, or both optimal and secondary habitats")

  if (is.null(species.dist$VALUE_1_prct)) warning(paste("Species", i, "has no secondary habitat in the study area"))
  if (is.null(species.dist$VALUE_2_prct)) warning(paste("Species", i, "has no primary habitat in the study area"))

  species.dist$Present <- 0
  species.dist[which(species.dist$VALUE_2_prct > threshold),"Present"] <- 1

  if (!opt.only) {
    species.dist[which(species.dist$VALUE_1_prct > threshold),"Present"] <- 1   # if using secondary habitats, add these as presences
  }

  return(species.dist[,c(1,ncol(species.dist))])
}


#' CONVERT SPECIES DATABASE FILES TO RASTERS
#' @param SPPPA `data.frame`  of species presense and absences
#'  with firt column being the cell ID ('PAGENAME') and the second the species
#'  presence/absence
#' @param mask.dir folder where study area mask rasters and dbf file with pixel IDs can be found
#'
#' @importFrom raster raster
#' @importFrom foreign read.dbf
#' @export
fun.dbf2raster <- function(SPPPA, mask.dir = NULL) {
  if (is.null(mask.dir)) stop("Must specify folder directory for mask files")

  if (any(!class(SPPPA) %in% "data.frame"))
    stop("Please make sure SPPPA is a data.frame")

  maskID <- read.dbf(list.files(path = mask.dir, full.names = TRUE, pattern = ".img.vat.dbf$"))
  maskk <- raster(x = list.files(path = mask.dir, full.names = TRUE, pattern = ".img$"))

  spp <- data.frame(PAGENAME = as.character(maskID$PageName), stringsAsFactors = FALSE)
  rownames(spp) <- spp$PAGENAME
  spp$val <- NA

  SPPPA[,1] <- as.character(SPPPA[,1])
  SPPPA[,2] <- as.numeric(as.character(SPPPA[,2]))
  rownames(SPPPA) <- SPPPA[,1]

  if (nrow(spp[SPPPA[,1],]) != length(SPPPA[,1])) stop("Cell IDs do not match")
  spp[SPPPA[,1],2] <- SPPPA[,2]

  xx <- maskk[]

  if (length(xx[!is.na(xx)]) != nrow(spp)) stop("Mask size inadequate")
  xx[!is.na(xx)] <- spp$val[xx[!is.na(xx)]]

  maskk[] <- xx
  return(maskk)
}


#' CONVERT LIST OF WEBS/SPP IN PIXELS INTO A PIXELS x SSP MATRIX
#' @param list.obj is a list of webs (if is.web = TRUE), or of vectors of species,
#'   per pixel that is to be converted into a matrix of pix X spp
#' @param is.web is `list.obj` a list of webs?
#'
#' @importFrom data.table rbindlist set dcast
#' @export
pixXspp_function <- function(list.obj = NULL, is.web = TRUE) {
  ## checks
  if (is.null(list.obj) |
      !is.list(list.obj) |
      !is.vector(list.obj))
    stop("Provide a list of webs or a vector of species")

  if (!is.list(list.obj) & is.web) {
    message("You provided a vector of species but 'is.web' is TRUE.\n
            'is.web' will be treated as FALSE")
    is.web <- FALSE
  }

  if (is.null(names(list.obj))) {
    stop("list.obj must be a named list.")
  }

  if (is.web) {
    ## make one column from aggregating all spp present in each pixel
    spp.ls <- data.table(Sp = unlist(apply(as.matrix(names(list.obj)), 1,
                                           function(x) rownames(list.obj[[x]]))))
    ## vector with number of spp per pixel
    no.spp <- apply(as.matrix(names(list.obj)), 1,
                    function(x) dim(list.obj[[x]])[1])
    ## associate page names by repeating pixel name as many times as there were spp.
    spp.ls[, PAGENAME := rep(names(list.obj), times = no.spp)]

    ## convert to presence/absence table, showing the number of times (max 1) a spp appears per pix
    master.scen <- with(spp.ls, table(PAGENAME, Sp))
    master.scen <- data.table(master.scen)

    ## check if everything's good
    errorPix <- master.scen[N > 1]$PAGENAME

    if (length(errorPix))
      stop("something's wrong. Pixels: ",
           paste(errorPix, collapse = ", "),
           " have some repeated species")

    ## convet to long format
    master.scen <- dcast(PAGENAME ~ Sp, data = master.scen, value.var = "N")

    ## remove column "1", that refers to empty pixels; add PAGENAME col
    master.scen$`1` <- NULL
  } else {
    ## make one column from aggregating all spp present in each pixel
    spp.ls <- lapply(list.obj, FUN = function(x) {
      if (is.vector(x) & !is.list(x)) {
        x <- data.table(Sp = x)
      } else {
        if (ncol(x) > 1) {
          stop("list.obj should be a list of vectors, or single columned tables")
        }
        names(x) = "Sp"
      }
      return(x)
    })
    spp.ls <- rbindlist(spp.ls, idcol = "PAGENAME")

    ## convert to table, showing the number of times (max 1) a spp appears per pix
    master.scen <- with(spp.ls, table(PAGENAME, Sp))
    master.scen <- as.data.frame.matrix(master.scen)
    master.scen <- data.table(master.scen, "PAGENAME" = rownames(master.scen))

    ## remove column "1", that refers to empty pixels; add PAGENAME col
    suppressWarnings(set(master.scen, NULL, j = "1", NULL))
  }
  return(master.scen)
}


#' CONVERT COLUMN NAMES IN MATRIX
#' @param M a `matrix` whose columns names will be changed
#' @param corresp is `data.frame` with rownames as the original column names
#'   and a single column with the final column names. If more than one original
#'   column matches a final one, the function will merge values using `merge.fun`
#' @param merge.fun funciton to use to merge columns. Default is `max`,
#'   which takes the maximum value across columns. Other possibilities are
#'   `mean` or `min`.
#' @param na.rm is passed to merge.fun
#' @export
col_convert <- function(M, corresp, merge.fun = "max", na.rm = TRUE) {
  ## subset habitats that have a correspondence and convert columns names
  M <- M[, which(colnames(M) %in% rownames(corresp))]    ## keep PAGENAME col
  colnames(M) <- corresp[colnames(M),]

  ## collapse habitats that are the same
  habs <- unique(colnames(M)[-1])
  M2 <- do.call(cbind.data.frame, lapply(habs, FUN = function(h) {
    if (dim(M[, colnames(M) %in% h, drop = FALSE])[2] > 1) {
      ## if there is more than one column of the same habitat, put NA/0 if both are NA/0, put 2 if one is 2, else put 1
      temp <- M[, colnames(M) %in% h, drop = FALSE]
      temp2 <- as.data.frame(apply(temp, 1, FUN = function(x) {
        if (all(is.na(x))) {return(NA)} else {
          return(
            eval(parse(text = paste0(
              merge.fun, "(x, na.rm = ", na.rm, ")"
            )))
          )
        }}))
      colnames(temp2) = h
    } else {
      temp2 <- M[, colnames(M) %in% h, drop = FALSE]
    }
    return(temp2)
  }))
}


#' EXTRACT PRESENCE ABSENCE MATRICES FROM SPECIES DISTRIBUTIONS
#'
#' @param files is a list of file paths to the distributions to be used
#' @param opt is either TRUE of FALSE (default) and determines if only optimal habitats are used as presences
#' @param k presence threshold \\[0,1\\], used if `opt = FALSE`. `opt` and
#'   `k` are only used for Luigi Maiorano's spp distributions
#' @export
do.PAmaster <- function(files, opt = FALSE, k = 0) {
  ## check file format
  if (!all(grepl(".img$", files)) & !all(grepl(".RData$", files))) {
    stop("\nDistribution files should all be .img or .RData files")
  }
  if (all(grepl(".img$", files))) {
    cat("\nConverting raster binary distributions to pixXspp matrix", append = TRUE)
    cat(paste0("\nStart:", date()), append = TRUE)
    ## Loading the spatial distribution of species
    spp.dist <- lapply(files, FUN = function(x) {
      if (grepl(".zip", folder)) {
        unzip(folder, files = x, exdir = sppdist.dir)
        temp <- raster(file.path(sppdist.dir, x))
      } else {
        temp <- raster(file.path(folder, x))
      }
      temp2 <- temp[]
      temp3 <- data.frame(X1 = rep(NA, length(temp2)), X2 = rep(NA, length(temp2)))
      temp3$X1[!is.na(mask10k[])] <- as.character(mask10kID$PageName)    # attributing same cellIDs that correspond to non-NAs in 10K grid
      temp3$X2 <- temp2
      colnames(temp3) = c("PAGENAME", toupper(sub("^m", "", sub("_.*", "", names(temp)))))
      temp3 <- temp3[!is.na(temp3$PAGENAME),]

      file.remove(file.path(sppdist.dir, x))

      return(temp3)
    })

    ## Make an pixel X spp matrix (at 10Km) by merging spp dataframes
    ## checks
    if (sum(!is.na(mask10k[])) != length(spp.dist[[1]]$PAGENAME)) {  # same number of cells?
      stop("number of pixels in SDMs differs for mask")
    }
    master <- Reduce(f = function(x, y) merge(x, y, all=TRUE, by = "PAGENAME", sort = FALSE), x = spp.dist)

    ## remove lakes and clean workspace
    master <- master[-lakes[!is.na(lakes[])],]

    cat("\nDone!", date(), append = TRUE)
    return(master)

  }
  if (all(grepl(".RData$", x))) {
    cat("\nConverting .RData distributions to pixXspp matrix", append = TRUE)
    cat(paste0("\nCalculating interaction matrix @ 10Km scale using a threshold higher than ", k, "% shared habitat"), append = TRUE)
    cat(paste0("\nStart:", date()), append = TRUE)

    ## Loading the spatial distribution of species
    spp.dist <- fun.load.many.rdata(file.list = x) # loads all spp spatial dist. in a list

    # replacing list names by spp names (removing "tab" prefix)
    names(spp.dist) <- gsub("(\\w)","\\U\\1", gsub("^tab_m", "", names(spp.dist)) , perl=TRUE)

    ## Make an pixel X spp matrix (at 10Km)
    master <- as.data.frame(matrix(data = NA,nrow = nrow(mask10kID), ncol = 1+length(spp.dist)))
    names(master) <- c("PAGENAME", names(spp.dist))
    master[,1] <- mask10kID$PageName

    for(i in 1:length(spp.dist)) {
      test <- fun.PRESENCE.ABSENCE(species.dist = spp.dist[[i]],threshold = k, opt.only = opt)
      cellID <- as.character(test$PAGENAME)

      master[master$PAGENAME %in% cellID,i+1] <- test[test$PAGENAME %in% cellID,"Present"]
    }

    cat("\nDone!", date(), append = TRUE)
    return(master)
  }
}

#' LOAD METRICS RESULTS FUNCTION
#' compiles simulation resuls and saves the data.table
#' @param bl.dir is the directory where baseline simulation results were saved
#' @template res.dir
#' @template out.dir
#' @param quant is the quantile threshold of extinction chosen
#' @param onlyBL load only baseline networks?
#' @param onlyScen load only scenario networks?
#'
#' @importFrom data.table data.table setkey rbindlist
#' @importFrom raster extension
#' @export
loadResultsMetrics <- function(bl.dir = NULL, res.dir = NULL, out.dir = NULL,
                               quant = c("No_ext_thresh", "min", "quant10", "quant25", "median", "quant50", "quant75", "quant90"),
                               onlyBL = FALSE, onlyScen = FALSE) {
  ## checks
  if (onlyBL & onlyScen) {
    stop("'onlyBL' and 'onlyScen' can't both be TRUE")
  }

  if ((onlyBL & is.null(bl.dir)) |
      (!onlyScen & is.null(bl.dir))) {
    stop("Please provide 'bl.dir'")
  }

  if ((onlyScen & is.null(res.dir)) |
      (!onlyBL & is.null(res.dir))) {
    stop("Please provide 'res.dir'")
  }

  quants <- c("No_ext_thresh", "min", "quant10", "quant25", "median", "quant50", "quant75", "quant90")
  if (!any(quant %in% quants)) {
    stop("Choose ONE of the following for 'quant':
         'No_ext_thresh', min', 'quant10', 'quant25', 'median', 'quant50', 'quant75', 'quant90'")
  }

  message("Compiling all_metrics")
  ## BASELINE NETWORKS
  folderName <- list.files(bl.dir, pattern = quant, full.names = TRUE)
  file.ls <- list.files(folderName, pattern = "metrics10k", full.names = TRUE)
  correctedFile.ls <- list.files(folderName, pattern = "metricsCorrect", full.names = TRUE)

  origFiles <- data.table(origFile = basename(file.ls), dir = as.character(dirname(file.ls)))
  corrFiles <- data.table(corrFile = basename(correctedFile.ls), dir = as.character(dirname(correctedFile.ls)))
  origFiles[, scen := sub(".*10kWdietFUND_", "", origFile)]
  corrFiles[, scen := sub(".*10kWdietFUND_", "", corrFile)]
  setkey(origFiles, dir, scen)
  setkey(corrFiles, dir, scen)
  file.lsDT <- corrFiles[origFiles]

  metrics_BL <- rbindlist(lapply(1:nrow(file.lsDT), FUN = function(x) {
    ## get the original and corrected files (if the latter exist compile.)
    origFile <- file.path(file.lsDT[x, dir], file.lsDT[x, origFile])
    corrFile <- file.path(file.lsDT[x, dir], file.lsDT[x, corrFile])
    if (!basename(corrFile) == "NA") {
      if (grepl("\\.R(d|D)ata", extension(basename(origFile)))) {
        metrics <- get(load(origFile))
      } else if (grepl("\\.rds", extension(basename(origFile)))) {
        metrics <- readRDS(origFile)
      } else stop("origFile needs to be .RData/.Rdata/.rds")

      if (grepl("\\.R(d|D)ata", extension(basename(corrFile)))) {
        metricsCorrected <- get(load(corrFile))
      } else if (grepl("\\.rds", extension(basename(corrFile)))) {
        metricsCorrected <- readRDS(corrFile)
      } else stop("corrFile needs to be .RData/.Rdata/.rds")

      ## change column names
      cols <- which(colnames(metrics) != "PAGENAME")
      colnames(metrics)[cols] <- paste0("BL_" , colnames(metrics)[cols])

      cols <- which(colnames(metricsCorrected) != "PAGENAME")
      colnames(metricsCorrected)[cols] <- sub("2", "", colnames(metricsCorrected)[cols]) %>%
        paste0("BL_" , .)

      ## remove old metrics and replace with corrected ones
      cols2Rm <- intersect(colnames(metrics)[2:length(metrics)], colnames(metricsCorrected)[2:length(metricsCorrected)])
      set(metrics, j = cols2Rm, value = NULL)

      setkey(metrics, PAGENAME)
      setkey(metricsCorrected, PAGENAME)
      metrics <- metrics[metricsCorrected]
    } else {
      if (grepl("\\.R(d|D)ata", extension(basename(origFile)))) {
        metrics <- get(load(origFile))
      } else if (grepl("\\.rds", extension(basename(origFile)))) {
        metrics <- readRDS(origFile)
      } else stop("origFile needs to be .RData/.Rdata/.rds")
      cols <- which(colnames(metrics) != "PAGENAME")
      colnames(metrics)[cols] <- paste0("BL_" , colnames(metrics)[cols])
    }

    metrics$SDM_stat <- if (grepl("*bin_GAM*", origFile)) "GAM" else "RF"
    return(metrics)
  }))
  metrics_BL$PAGENAME <- as.character(metrics_BL$PAGENAME)

  if (onlyBL) {
    message("Returning metrics from baseline networks only")
    return(metrics_BL)
  }

  ## SCENARIO NETWORKS
  ## recursive list.files was stalling
  ## needs to come first for memory reasons
  file.ls <- do.call(c, lapply(list.files(res.dir, full.names = TRUE), function(x) {
    do.call(c, lapply(list.files(x, pattern = quant, full.names = TRUE), function(xx) {
      list.files(xx, pattern = "metrics10k", full.names = TRUE)
    }))
  }))

  correctedFile.ls <- do.call(c, lapply(list.files(res.dir, full.names = TRUE), function(x) {
    do.call(c, lapply(list.files(x, pattern = quant, full.names = TRUE), function(xx) {
      list.files(xx, pattern = "metricsCorrect", full.names = TRUE)
    }))
  }))

  origFiles <- data.table(origFile = basename(file.ls), dir = as.character(dirname(file.ls)))
  corrFiles <- data.table(corrFile = basename(correctedFile.ls), dir = as.character(dirname(correctedFile.ls)))
  origFiles[, scen := sub(".*10kWdietFUND_", "", origFile)]
  corrFiles[, scen := sub(".*10kWdietFUND_", "", corrFile)]
  setkey(origFiles, dir, scen)
  setkey(corrFiles, dir, scen)
  file.lsDT <- corrFiles[origFiles]

  all_metrics <- rbindlist(lapply(1:nrow(file.lsDT), FUN = function(x) {
    ## get the original and corrected files (if the latter exist compile.)
    origFile <- file.path(file.lsDT[x, dir], file.lsDT[x, origFile])
    corrFile <- file.path(file.lsDT[x, dir], file.lsDT[x, corrFile])
    if (!basename(corrFile) == "NA") {
      if (grepl("\\.R(d|D)ata", extension(basename(origFile)))) {
        metrics <- get(load(origFile))
      } else if (grepl("\\.rds", extension(basename(origFile)))) {
        metrics <- readRDS(origFile)
      } else stop("origFile needs to be .RData/.Rdata/.rds")

      if (grepl("\\.R(d|D)ata", extension(basename(corrFile)))) {
        metricsCorrected <- get(load(corrFile))
      } else if (grepl("\\.rds", extension(basename(corrFile)))) {
        metricsCorrected <- readRDS(corrFile)
      } else stop("corrFile needs to be .RData/.Rdata/.rds")

      ## change column names
      colnames(metricsCorrected)[2:length(metricsCorrected)] <- sub("2", "", colnames(metricsCorrected)[2:length(metricsCorrected)])

      ## remove old metrics and replace with corrected ones
      cols2Rm <- intersect(colnames(metrics)[2:length(metrics)], colnames(metricsCorrected)[2:length(metricsCorrected)])
      set(metrics, j = cols2Rm, value = NULL)

      setkey(metrics, PAGENAME)
      setkey(metricsCorrected, PAGENAME)
      metrics <- metrics[metricsCorrected]

    } else {
      if (grepl("\\.R(d|D)ata", extension(basename(origFile)))) {
        metrics <- get(load(origFile))
      } else if (grepl("\\.rds", extension(basename(origFile)))) {
        metrics <- readRDS(origFile)
      } else stop("origFile needs to be .RData/.Rdata/.rds")
    }
    ## remove general scenario simulation parent directory
    origFile <- sub(res.dir, "", origFile, fixed = TRUE)

    metrics$SDM_stat <- if (grepl("*bin_GAM*", origFile)) "GAM" else "RF"
    metrics$LUC <- if (grepl("LUC_", origFile) & !grepl("noLUC", origFile))
      sub(".*SSP", "SSP", sub("_allhab.*(Rdata|RData|rds)$", "", basename(origFile))) else "noLUC"
    metrics$GCM <- if (grepl("CC_", origFile))
      sub(".*dietFUND_", "", sub("_rcp.*", "", basename(origFile))) else "noCC"
    metrics$RCP <- if (grepl("CC_", origFile))
      sub(".*dietFUND_.*_rcp", "rcp", sub("_wm_bin.*", "", basename(origFile))) else "noCC"
    metrics$Invasions <- if (grepl("noInvs", origFile) |
                             (grepl("noCC", origFile) | grepl("current", origFile)))
      "noInvs" else "Invs"
    metrics$IUCN <- if (!grepl("noIUCN", origFile))
      sub("(.*)(CR|CR_EN|CR_EN_VU)", "\\2", sub("_allhab.*(Rdata|RData|rds)$", "", basename(origFile))) else "noIUCNext"
    metrics$resDir <- file.lsDT[x, dir]

    return(metrics)
  }))
  all_metrics$PAGENAME <- as.character(all_metrics$PAGENAME)

  if (onlyScen) {
    message("Returning metrics from scenario networks only - careful, some maybe be 'new' networks due to SDMs" )
    return(all_metrics)
  }

  setkey(all_metrics, PAGENAME, SDM_stat)
  setkey(metrics_BL, PAGENAME, SDM_stat)
  all_metrics <- metrics_BL[all_metrics]

  ## REPLACE FUTURE NAs  in S by 0s.
  ## some future network metrics have NA's in S because no metrics were calculated
  ## even though Pext and Sext were. These networks should've had S = 0 (and other metrics)
  all_metrics[!is.na(BL_S) & is.na(S), `:=` (S = 0,
                                             L = 0,
                                             C = 0,
                                             propB = 0,
                                             propI = 0,
                                             propT = 0,
                                             mean.TL = 0,
                                             Invasiv = 0)]
  all_metrics[!is.na(BL_S) & S == 0 & is.na(Pext), Pext := BL_S]

  ## Creating robustness variable
  all_metrics$Robust <- all_metrics$Sext/all_metrics$BL_No.SecCons

  message("Saving all_metrics...")
  saveRDS(all_metrics, file = file.path(out.dir, paste0("all_metrics", quant, ".rds")))
  message("Done!")

  return(all_metrics)
}


#' LOAD BASELINE MASTER MATRICES (SPP X PIXEL MATRICES) FUNCTION
#'
#' compiles simulation results and saves the data.table
#'
#' @param file.ls is the list of file.paths for the master matrices
#' @template res.dir
#' @template out.dir
#' @param quant is the quantile threshold of extinction chosen
#' @param useCache is NULL, but will default to TRUE if the argument is not defined in the parent.frame()
#' @param cacheRepo passed to `reproducible::Cache`.
#' @template dietcat
#'
#' @importFrom reproducible Cache
#' @export
compileResultsMasterBL <- function(file.ls, res.dir, out.dir, quant, dietcat,
                                   useCache = NULL, cacheRepo = options("reproducible.cachePath")) {
  if (is.null(useCache)) {
    useCache <- tryCatch(get("useCache", envir = parent.frame()),
                         error = function(e) TRUE)
  }

  if (is.null(cacheRepo)) {
    cacheRepo <- tryCatch(get("cacheRepo", envir = parent.frame()),
                          error = function(e) tempdir())
  }

  message("Compiling all_masterBL")

  all_masterBL <- Cache(.compileMasterBL,
                        file.ls = file.ls,
                        dietcat = dietcat,
                        useCache = useCache,
                        cacheRepo = cacheRepo,
                        userTags = c("compileResultsMaster", "all_masterBL"),
                        omitArgs = c("userTags"))

  ## make PAGENAME a character
  all_masterBL[, PAGENAME := as.character(PAGENAME)]

  ## replace NAs with 0s (absences)
  all_masterBL <- replaceNAs.data.table(all_masterBL, )

  if (!file.exists(file.path(out.dir, paste0("all_masterBL", quant, ".rds"))) |
      !useCache | useCache == "overwrite") {
    suppressMessages(dir.create(out.dir, recursive = TRUE))
    message("Saving all_masterBL...")
    saveRDS(all_masterBL, file = file.path(out.dir, paste0("all_masterBL", quant, ".rds")))
    message("Done!")
  }

  ## clean ws
  return(all_masterBL)
}


#' LOAD SCENARIO MASTER MATRICES (SPP X PIXEL MATRICES) FUNCTION
#' compiles simulation results and saves the data.table
#'
#' @param file.ls is the list of file.paths for the master matrices
#' @template res.dir
#' @template out.dir
#' @param quant is the quantile threshold of extinction chosen
#' @param useCache is NULL, but will default to TRUE if the argument is not defined in the parent.frame()
#' @param cacheRepo is NULL, but will default to temp.dir() if the argument is not defined in the parent.frame()
#' @template dietcat
#'
#' @importFrom reproducible Cache
#' @export
compileResultsMasterScen <- function (file.ls, res.dir, out.dir, quant, dietcat,
                                      useCache = NULL, cacheRepo = NULL) {
  if (is.null(useCache)) {
    useCache <- tryCatch(get("useCache", envir = parent.frame()),
                         error = function(e) TRUE)
  }

  if (is.null(cacheRepo)) {
    cacheRepo <- tryCatch(get("cacheRepo", envir = parent.frame()),
                          error = function(e) tempdir())
  }

  ## SCENARIO NETWORKS
  ## recursive list.files was stalling
  ## needs to come first for memory reasons
  message("Compiling all_masterScen")

  all_masterScen <- Cache(.compileMasterScen,
                          file.ls = file.ls,
                          dietcat = dietcat,
                          useCache = useCache,
                          cacheRepo = cacheRepo,
                          userTags = c("compileResultsMaster", "all_masterScen"),
                          omitArgs = c("userTags"))

  ## make PAGENAME a character
  all_masterScen[, PAGENAME := as.character(PAGENAME)]

  ## replace NAs with 0s (absences)
  all_masterScen <- replaceNAs.data.table(all_masterScen)


  out.dir <- file.path(res.dir, "Analyses")
  if (!file.exists(file.path(out.dir, paste0("all_masterScen", quant, ".rds"))) |
      !useCache | useCache == "overwrite") {
    suppressMessages(dir.create(out.dir, recursive = TRUE))
    message("Saving all_masterScen...")
    saveRDS(all_masterScen, file = file.path(out.dir, paste0("all_masterScen", quant, ".rds")))
    message("Done!")
  }

  return(all_masterScen)
}


#' internal function
#'
#' @param file.ls file list
#' @template dietcat
#'
#' @importFrom data.table rbindlist
.compileMasterBL <- function(file.ls, dietcat) {
  rbindlist(fill = TRUE,
            use.names = TRUE,
            l = lapply(X = file.ls,
                       FUN = makeAndSaveMasterBL,
                       redo = FALSE,
                       returnMaster = TRUE,
                       save = FALSE))
}

#' internal function
#'
#' @param file.ls file list
#' @template dietcat
#'
#' @importFrom data.table rbindlist
.compileMasterScen <- function(file.ls, dietcat) {
  rbindlist(fill = TRUE,
            use.names = TRUE,
            l = lapply(X = file.ls,
                       FUN = makeAndSaveMasterScen,
                       redo = FALSE,
                       returnMaster = TRUE,
                       save = FALSE))
}


#' COMPILE SDM PROJECTIONS PRESENCE/ABSENCES
#'
#' @param fileLs is the list of file names (not paths!) to the different SDM projections
#' @param fileFolder is the path to the parent folder containing all files (and directories to the SDM projections)
#'
#' @importFrom data.table rbindlist dcast as.data.table
#' @export
compileSDMsppRichness <- function(fileLs, fileFolder) {
  outputList <- lapply(fileLs,
                       FUN = function(f){
                         load(file.path(fileFolder, f))
                         master <- get(grep("^master\\.", ls(), value = TRUE)) %>%
                           as.data.table(.)
                         ## make sure lakes are removed (shouldn't change the data)
                         lakePix <- mask10kID$PageName[mask10k[!is.na(lakes)]]
                         master <- master[!PAGENAME %in% lakePix]
                         ## calculate original S
                         cols <- grep("PAGENAME", names(master), invert = TRUE)
                         master[, S := rowSums(.SD), .SDcols = cols]

                         ## make a scenario ID column
                         master[, CCscen := sub("pixXspp_10k_(.*)_wm.*", "\\1", sub(".Rdata", "", f))]
                         master[, SDMstat := sub(".*bin_(.*)_all.*", "\\1", sub(".Rdata", "", f))]

                         return(master[, .(PAGENAME, S, SDMstat, CCscen)])
                       })

  sppRichness <- rbindlist(outputList, use.names = TRUE, fill = TRUE)

  sppRichness <- dcast(sppRichness, SDMstat + PAGENAME ~ CCscen, value.var = "S") %>%
    melt(., id.vars = c("SDMstat", "PAGENAME", "current"), variable.name = "CCscen", value.name = "Sfuture")

  sppRichness <- as.data.table(sppRichness)
  return(sppRichness)
}



#' MAKE AND SAVE BASELINE MASTER MATRICES (SPP X PIXEL MATRICES) FUNCTION
#' compiles all scenario networks into a data.table and saves it to an
#'  .rds file
#' @param x is the .RData file containing the list of networks.
#' @param returnMaster controls whether the function returns the master matrix or not. Defaults
#'   to FALSE and returns the saved file path
#' @param redo controls whether masters will be recalculated and re-saved
#'   in case their saved files already exist
#' @param save controls saving
#' @param dietcat diet categories (i.e., lowest trophic level nodes that are ubiquitous across all networks)
#'
#' @importFrom raster extension
#' @export
makeAndSaveMasterBL <- function(x, returnMaster = FALSE, dietcat,
                                redo = FALSE, save = TRUE) {
  ## checks
  if (!save & returnMaster)
    warning("neither saving nor returning master are activated.
            The function will return the master matrix")

  if (grepl("\\.R(D|d)ata", extension(basename(x)))) {
    outFile <- sub("\\.R(D|d)ata", "\\.rds", sub("spp10kW", "master", basename(x)))
  } else if (grepl("\\.rds", extension(basename(x)))) {
    outFile <- sub("spp10kW", "master", basename(x))
  } else {
    stop(paste(x, "needs to be .RData/.Rdata/.rds"))
  }
  outFile <- file.path(dirname(x), outFile)

  if (file.exists(outFile) & isFALSE(redo)) {
    message(paste(outFile, "already exists and 'redo' is", redo,
                  "\nskipping..."))
    if (returnMaster) {
      master <- readRDS(outFile)
      return(master)
    } else
      return(outFile)
  } else {
    if (grepl("\\.R(D|d)ata", extension(basename(x)))) {
      load(x)
    } else {
      pixelXspp.ls <- readRDS(x)
    }

    master <- pixXspp_function(list.obj = pixelXspp.ls)
    rm(pixelXspp.ls); for (i in 1:10) gc()
    ## remove diet categories
    keepCols <- setdiff(names(master), dietcat)
    master <- master[, ..keepCols]

    master$SDM_stat <- if (grepl("*bin_GAM*", x)) "GAM" else "RF"

    ## save
    if (save) {
      saveRDS(master, outFile)
      if (returnMaster)
        return(master) else
          return(outFile)
    } else {
      return(master)
    }
  }
}


#' MAKE AND SAVE SCENARIO MASTER MATRICES (SPP X PIXEL MATRICES) FUNCTION
#' compiles all scenario networks into a sparse matrix and saves it to an
#'  .rds file
#' @param x is the .RData file containing the list of networks.
#' @param returnMaster controls whether the function returns the master matrix or not. Defaults
#'   to FALSE and returns the saved file path
#' @param redo controls whether masters will be recalculated and re-saved
#'   in case their saved files already exist
#' @param save controls saving
#' @template dietcat
#'
#' @importFrom raster extension
#' @export
makeAndSaveMasterScen <- function(x, returnMaster = FALSE, dietcat,
                                  redo = FALSE, save = TRUE) {
  ## checks
  if (!save & returnMaster)
    warning("neither saving nor returning master are activated.
            The function will return the master matrix")

  if (grepl("\\.R(D|d)ata", extension(basename(x)))) {
    outFile <- sub("\\.R(D|d)ata", "\\.rds", sub("spp10kW", "master", basename(x)))
  } else if (grepl("\\.rds", extension(basename(x)))) {
    outFile <- sub("spp10kW", "master", basename(x))
  } else {
    stop(paste(x, "needs to be .RData/.Rdata/.rds"))
  }

  outFile <- file.path(dirname(x), outFile)

  if (file.exists(outFile) & isFALSE(redo)) {
    message(paste(outFile, "already exists and 'redo' is", redo,
                  "\nskipping..."))

    if (returnMaster) {
      master <- readRDS(outFile)
      return(master)
    } else
      return(outFile)
  } else {
    if (grepl("\\.R(D|d)ata", extension(basename(x)))) {
      load(x)
    } else {
      pixelXspp.ls <- readRDS(x)
    }
    master <- pixXspp_function(list.obj = pixelXspp.ls)
    rm(pixelXspp.ls); for (i in 1:10) gc()

    ## remove diet categories
    keepCols <- setdiff(names(master), dietcat)
    master <- master[, ..keepCols]

    ## keep only filename
    xx <- basename(x)

    master$SDM_stat <- if (grepl("*bin_GAM*", xx)) "GAM" else "RF"
    master$LUC <- if (grepl("LUC_", xx) & !grepl("noLUC", xx))
      sub("(.*_)(SSP[[:digit:]]{1}_[[:digit:]]{4})(_.*)", "\\2", xx) else "noLUC"
    master$GCM <- if (grepl("CC_", xx))
      sub(".*dietFUND_", "", sub("_rcp.*", "", xx)) else "noCC"
    master$RCP <- if (grepl("CC_", xx))
      sub("(.*_)(rcp[[:digit:]]{1,2})(_.*)", "\\2", xx) else "noCC"
    master$Invasions <- if (grepl("noInvs", xx) |
                            (grepl("noCC", xx) | grepl("current", xx)))
      "noInvs" else "Invs"
    master$IUCN <- if (!grepl("noIUCN", xx))
      sub("(.*)(CR|CR_EN|CR_EN_VU)", "\\2", sub("_allhab.*(Rdata|RData|rds)$", "", xx)) else "noIUCNext"

    ## save
    if (save) {
      saveRDS(master, outFile)
      if (returnMaster)
        return(master) else
          return(outFile)
    } else {
      return(master)
    }
  }
}


#' LOAD AND SUMMARISE METRICS RESULTS FUNCTION
#'
#' compiles and summarises simulation network metrics results.
#'   this function is a wrapper to `loadResultsMetrics` and `calcNetworkMetricsSummStats`
#'   functions and passes all arguments to these functions, as well as `reproducible::Cache`
#'
#' @param ... passed to `loadResultsMetrics`, `calcNetworkMetricsSummStats`
#'   and `reproducible::Cache`
#'
#' @importFrom reproducible Cache
#' @export
loadAndSummarizeResults <- function(...) {
  dots <- list(...,
               userTags = get("userTags", envir = parent.frame()),
               cacheRepo = get("cacheRepo", envir = parent.frame()))

  if (is.null(dots$cacheRepo))
    dots$cacheRepo <- getOption("reproducible.cachePath")
  all_metrics <- Cache(loadResultsMetrics, bl.dir = dots$bl.dir,
                       res.dir = dots$res.dir, out.dir = dots$out.dir,
                       quant = dots$quant,
                       cacheRepo = dots$cacheRepo, useCache = dots$useCache,
                       userTags = dots$userTags)

  ## transform Robustness
  all_metrics[, invRobust := 1-Robust]

  ## calculate link density:
  all_metrics[, BL_LD := BL_L/BL_S]
  all_metrics[, LD := L/S]

  ## calculate deltas
  all_metrics[, deltaS := S - BL_S]
  all_metrics[, deltaL := L - BL_L]
  all_metrics[, deltaLD := LD - BL_LD]
  all_metrics[, deltaC := C - BL_C]
  all_metrics[, absDeltaC := abs(C - BL_C)]
  all_metrics[, deltaQ := Q - BL_Q]
  all_metrics[, deltaSDnormGen := SDnormGen - BL_SDnormGen]
  all_metrics[, deltaSDnormVul := SDnormVul - BL_SDnormVul]
  all_metrics[, deltapropB := propB - BL_propB]
  all_metrics[, deltapropI := propI - BL_propI]
  all_metrics[, deltapropT := propT - BL_propT]
  all_metrics[, deltapropOmn := propOmn - BL_propOmn]
  all_metrics[, deltamean.TL := mean.TL - BL_mean.TL]
  all_metrics[, deltamax.TL := max.TL - BL_max.TL]
  all_metrics[, deltasd.TL := sd.TL - BL_sd.TL]

  # if (is.null(dots$measure.vars)) {
  #   dots$measure.vars <- grep("PAGENAME|SDM_stat|LUC$|PA$|GCM|RCP|Invasions|resDir|code",
  #                             names(all_metrics), invert = TRUE, value = TRUE)
  # }

  all_metrics <- calcNetworkMetricsSummStats(DT = all_metrics,
                                             stats = dots$stats,
                                             id.vars = dots$id.vars,
                                             measure.vars = dots$measure.vars)
  all_metrics
}


#' LOAD TEMPORAL BETA DIVERSITY RESULTS FUNCTION
#'
#' compiles calculations of temporal beta diversity results and outputs the data.table
#'
#' @param file.ls list of files to compile
#'
#' @importFrom data.table rbindlist
#' @export
loadResultsTempBetaDiv <- function(file.ls) {
  message("Compiling temporal beta diversity")
  Tbeta_div <- rbindlist(use.names = TRUE,
                         l = lapply(file.ls, FUN = function(x) {
                           tempTaxonBetaDiv <- readRDS(x)
                           ## remove general scenario simulation parent directory
                           basedir <- sub("Analyses/.*", "Analyses/", dirname(x))
                           origFile <- sub(basedir, "", x, fixed = TRUE)

                           tempTaxonBetaDiv$Taxon <- if (!grepl("mammal|reptile|bird|amphibian", origFile))
                             "All" else
                               sub(".*(mammal|reptile|bird|amphibian).*", "\\1", origFile)

                           tempTaxonBetaDiv$SDM_stat <- if (grepl("*_GAM_*", origFile)) "GAM" else "RF"
                           tempTaxonBetaDiv$LUC <- if (grepl("noLUC", origFile) & !grepl("SSP", origFile))
                             "noLUC" else sub("(.*_)(SSP[[:digit:]]{1}_[[:digit:]]{4})(_.*)", "\\2", origFile)
                           tempTaxonBetaDiv$GCM <- if (grepl("noCC", origFile) | grepl("current", origFile))
                             "noCC" else sub(".*BetaDiv_(mammal_|reptile_|bird_|amphibian_)?", "", sub("_rcp.*", "", origFile))
                           tempTaxonBetaDiv$RCP <- if (grepl("noCC", origFile) | grepl("current", origFile))
                             "noCC" else sub("(.*_)(rcp[[:digit:]]{1,2})(_.*)", "\\2", origFile)
                           tempTaxonBetaDiv$Invasions <- if (grepl("noInvs", origFile) |
                                                             (grepl("noCC", origFile) | grepl("current", origFile)))
                             "noInvs" else "Invs"
                           tempTaxonBetaDiv$IUCN <- if (grepl("noInvs", origFile) |
                                                        (grepl("noCC", origFile) | grepl("current", origFile)))
                             sub("(.*)(CR|CR_EN|CR_EN_VU)", "\\2", sub("_allhab.*(Rdata|RData|rds)$", "", origFile)) else "noIUCNext"
                           tempTaxonBetaDiv$resDir <- sub("/TaxonBeta.*", "", origFile)

                           return(tempTaxonBetaDiv)
                         }))

  cols <- grep("beta|diss", names(Tbeta_div), value = TRUE)
  Tbeta_div[, lapply(.SD, as.numeric), .SDcols = cols]
  return(Tbeta_div)
}


#' LOAD SPATIAL BETA DIVERSITY RESULTS FUNCTION
#'
#' compiles calculations of spatial beta diversity results and outputs the data.table
#'
#' @inheritParams loadResultsTempBetaDiv
#'
#' @importFrom data.table rbindlist
#' @export
loadResultsSpaceBetaDiv <- function(file.ls) {
  message("Compiling spatial beta diversity")
  Sbeta_div <- rbindlist(use.names = TRUE,
                         l = lapply(file.ls, FUN = function(x) {
                           spaceTaxonBetaDiv <- readRDS(x)

                           if (is(spaceTaxonBetaDiv, "list")) {
                             if (all(sapply(spaceTaxonBetaDiv, function(x) is.list(x)))) {
                               spaceTaxonBetaDivMean <- sapply(spaceTaxonBetaDiv, function(x) {
                                 betaID <- grep("(B|b)eta", names(x))
                                 data.table(meanBetadiversity = x[[betaID]])
                               }, simplify = FALSE)

                               spaceTaxonBetaDivMean <- rbindlist(spaceTaxonBetaDivMean, use.names = TRUE, idcol = "ID")

                             } else if (all(sapply(spaceTaxonBetaDiv, function(x) is(x, "dist")))) {
                             spaceTaxonBetaDivMean <- rbindlist(lapply(spaceTaxonBetaDiv, FUN = calcNetMeanDistances),
                                                                idcol = "ID")
                             spaceTaxonBetaDivMean[, rep := as.numeric(as.factor(rep))]  ## make sure its a number
                               } else stop("Expecting a list of 'dist' objects or of lists of",
                                           " alpha, beta and gamma diversity values")
                           } else if (is(spaceTaxonBetaDiv, "dist")) {
                             spaceTaxonBetaDivMean <- calcNetMeanDistances(spaceTaxonBetaDiv)
                             spaceTaxonBetaDivMean[, rep := 1]
                           } else stop("Expecting a list (of 'dist' objects or of lists of",
                                       " alpha, beta and gamma diversity values) or",
                                       " a 'dist' object.")

                           ## remove general scenario simulation parent directory
                           basedir <- sub("Analyses/.*", "Analyses/", dirname(x))
                           origFile <- sub(basedir, "", x, fixed = TRUE)

                           spaceTaxonBetaDivMean$SDM_stat <- if (grepl("*_GAM_*", origFile)) "GAM" else "RF"
                           spaceTaxonBetaDivMean$LUC <- if (grepl("noLUC", origFile) & !grepl("SSP", origFile))
                             "noLUC" else sub("(.*_)(SSP[[:digit:]]{1}_[[:digit:]]{4})(_.*)", "\\2", origFile)
                           spaceTaxonBetaDivMean$GCM <- if (grepl("noCC", origFile) | grepl("current", origFile))
                             "noCC" else sub(".*BetaDiv_(mammal_|reptile_|bird_|amphibian_)?", "", sub("_rcp.*", "", origFile))
                           spaceTaxonBetaDivMean$RCP <- if (grepl("noCC", origFile) | grepl("current", origFile))
                             "noCC" else sub("(.*_)(rcp[[:digit:]]{1,2})(_.*)", "\\2", origFile)
                           spaceTaxonBetaDivMean$Invasions <- if (grepl("noInvs", origFile) |
                                                                  (grepl("noCC", origFile) | grepl("current", origFile)))
                             "noInvs" else "Invs"
                           spaceTaxonBetaDivMean$IUCN <- if (!grepl("noIUCN", origFile))
                             sub("(.*)(CR|CR_EN|CR_EN_VU)", "\\2", sub("_allhab.*(Rdata|RData|rds)$", "", origFile)) else "noIUCNext"
                           spaceTaxonBetaDivMean$resDir <- sub("/Taxon.*Beta.*", "", origFile)

                           return(spaceTaxonBetaDivMean)
                         }))

  return(Sbeta_div)
}
