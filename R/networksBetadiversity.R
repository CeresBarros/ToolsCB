#' ----------------------------------------------
#' NETWORKS BETADIVERSITY
#' ----------------------------------------------

#' TEMPORAL BETA-DIVERSITY PER TETRAPOD FAMILY WRAPPER FUNCTION   ---------------------------
#'
#' This function loads and prepares networks to calculate temporal
#' spp turnover (beta diversity) for each tetrapod family separately and saves the outputs.
#' beta diversity is calculated per pairs of pixels of the same location,
#' but different scenarios (a baseline/current scenario and a change scenario).
#' this function was designed to be parallelelized.

#' @param ff is the a file path to list of networks obtained from a scenario of change.
#' @param BLweb a baseline web (not a file path) to be compared with the scenario webs
#' @param res.dir is the "root" directory where "scenario networks" were stored and out.dir is
#'   were results will be placed (in a folder tree respecting scenario folders.)
#' @template dietcat
#' @param families determines for which families beta div will be calculated (this subsets the species)
#' @param toDo controls if beta-diversity has to be recalculated for 'all' networks or
#'   just 'missing' ones.
#' @param method is passed to \code{networkTempBetaDiv}
#' @param mode is passed to \code{networkTempBetaDiv}
#' @param cacheRepo passed to \code{reproducible::Cache}
#'
#' @importFrom crayon red
#' @importFrom tools file_path_sans_ext
#' @importFrom reproducible Cache
#' @importFrom data.table rbindlist
#'
#' @export
calcTempBetaDivBAMR <- function(ff, BLweb, res.dir, out.dir, dietcat = dietcat,
                                families = c("allTaxon", "bird", "mammal", "reptile", "amphibian"),
                                toDo = "all", method = "all", mode = "composition",
                                cacheRepo = options("reproducible.cachePath")) {
  ## checks
  if (any(!families %in% c("allTaxon", "bird", "mammal", "reptile", "amphibian")))
    stop("'families' must be one, or several of
         c('allTaxon', 'bird', 'mammal', 'reptile', 'amphibian')")

  if (any(mode %in% c("decompHills", "pairwiseHills")) &
      any(families != "allTaxon")) {
    message(red(paste("Subsetting networks by species family will probably",
                      "result in diconnected nodes, which are problematic when",
                      "calculating link and species turnover with mode", mode)))
  }

  ## change "allTaxon" to nothing for file suffixes and subset
  ## selected taxon
  familiesSufs <- c("", "bird", "mammal", "reptile", "amphibian")
  names(familiesSufs) <- c("allTaxon", "bird", "mammal", "reptile", "amphibian")
  familiesSufs <- familiesSufs[families]

  cat(paste0("***\n", "Start\n", date(), ":\n", ff, "\n"), append = TRUE)

  ## make output file suffix and directory
  scen <- sub(res.dir, "", dirname(ff))
  fileSuf <- file_path_sans_ext(basename(ff))
  fileSuf <- sub("spp10kWdietFUND_", "", fileSuf)

  out.folder <- file.path(out.dir, scen)
  if (!dir.exists(out.folder))
    dir.create(out.folder, recursive = TRUE)

  ## make all file names and check for their existence
  ## only load webs if some files need to be created, or all computations need to be repeated
  outFiles <- if (mode == "composition") {
    # file.path(out.folder, paste0("TaxonTempBetaDiv_", fileSuf, ".rds"))
    file.path(out.folder, paste0("TaxonTempBetaDiv_", familiesSufs, "_", fileSuf, ".rds"))
  } else {
    # file.path(out.folder, paste0("SppLinkTempBetaDiv_", fileSuf, ".rds"))
    file.path(out.folder, paste0("SppLinkTempBetaDiv_", familiesSufs, "_", fileSuf, ".rds"))
  }

  outFiles <- sub("__", "_", outFiles) ## remove double "__" when allTaxon is present
  names(outFiles) <- families
  rm(familiesSufs) ## clean-up to avoid confusion

  ## check which files need to be done, skip files done by subsetting families.
  if (any(file.exists(outFiles)) & toDo == "missing") {
    cat(paste0("***\n", "Already done:\n",
               paste(outFiles[file.exists(outFiles)], collapse = "\n"),
               "\n... skipping\n"), append = TRUE)
    outFiles <- outFiles[!file.exists(outFiles)]
    families <- names(outFiles)
  }

  if (length(families)) {
    for (fam in families) {
      outFile <- outFiles[fam]
      cat(paste0("***\n", "Calculating:\n", outFile))

      charMatch <- if (fam != "allTaxon") {
        paste0("^",  toupper(substr(fam, 1, 1)), "[[:xdigit:]]{1}")
      } else NULL

      if (grepl("\\.R(D|d)ata", extension(basename(ff)))) {
        tryCatch(load(ff), error = function(e) NULL)
      } else if (grepl("\\.rds", extension(basename(ff)))) {
        tryCatch(pixelXspp.ls <- readRDS(ff), error = function(e) NULL)
      } else stop("ff must be a .RData/.Rdata/.rds file")

      if (!exists("pixelXspp.ls")) {
        cat(paste0("***\n", "Loading failed for:\n", ff, "\n"), append = TRUE)
      } else {
        ## TEMPORAL BETA_DIVERSITY CALCULATION ---------------------------
        ## iterate through pixels (i.e, networks) and subset networks to the family's species
        tempBetaDiv <- Cache(Map,
                             f = .networkTempBetaDiversityFamilies,
                             i = names(pixelXspp.ls),
                             MoreArgs = list(web1.ls = BLweb,
                                             web2.ls = pixelXspp.ls,
                                             dietcat = dietcat,
                                             charMatch = charMatch,
                                             method = method,
                                             mode = mode),
                             .cacheExtra = c(ff, length(pixelXspp.ls),
                                             object.size(pixelXspp.ls),
                                             object.size(BLweb),
                                             charMatch, method, mode),
                             cacheRepo = cacheRepo,
                             userTags = c("tempTaxonBetaDiv"),
                             omitArgs = c("userTags", "MoreArgs"))
        tempBetaDiv <- rbindlist(tempBetaDiv, use.names = TRUE, fill = TRUE)

        ## save data.table
        cat("***\n", "Saving...\n", append = TRUE)
        saveRDS(tempBetaDiv, file = outFile)
        cat(paste0("Saved!\n", date(), "\n"), append = TRUE)

        return(outFile)
      }
    }
  }
}


#' NETWORK BETA DIVERSITY PER FAMILY USING WEBS AS INPUTS
#'
#' Internal function used to calculate pairwise beta-diversity
#' for two networks per species family. The network pairs to compare should have the same
#' positions in the network lists supplied - e.g. the first network in list 1 is
#' compared with the first network in list 2, the second network in list 1 is
#' compared with the second network in list 2, and so on.
#'
#' @param i is the index to subset network lists by (e.g. pixel name if calculating temporal beta-diversity)
#' @param web1.ls is a list of networks to be indexed by i
#' @param web2.ls is a list of networks to be indexed by i
#' @template dietcat
#' @param charMatch is used to subset species in the networks beloging to a specific family.
#' @param method is passed to \code{networkTempBetaDiv}
#' @param mode is passed to \code{networkTempBetaDiv}

.networkTempBetaDiversityFamilies <- function(i, web1.ls, web2.ls, dietcat,
                                              charMatch, method, mode = "composition") {
  # if (i == "AU366") browser()
  # print(i)
  ## extract networks of pixel i from lists
  web1 <- as.matrix(web1.ls[[i]])
  web2 <- as.matrix(web2.ls[[i]])

  ## get row and column names in networks that correspond to the species
  ## in the current family
  ## keep diet categories -- necessary to calculate link turnover
  if (!is.null(charMatch)) {
    rowsWeb1 <- grep(paste0(charMatch, "|", paste(dietcat, collapse = "|")), rownames(web1), value = TRUE)
    colsWeb1 <- grep(paste0(charMatch, "|", paste(dietcat, collapse = "|")), colnames(web1), value = TRUE)

    rowsWeb2 <- grep(paste0(charMatch, "|", paste(dietcat, collapse = "|")), rownames(web2), value = TRUE)
    colsWeb2 <- grep(paste0(charMatch, "|", paste(dietcat, collapse = "|")), colnames(web2), value = TRUE)
  } else {
    rowsWeb1 <- rownames(web1)
    colsWeb1 <- colnames(web1)

    rowsWeb2 <- rownames(web2)
    colsWeb2 <- colnames(web2)
  }

  ## subset webs
  web1 <- web1[rowsWeb1, colsWeb1, drop = FALSE]
  web2 <- web2[rowsWeb2, colsWeb2, drop = FALSE]

  ## add missing diet categories  with 0s.
  webLs <- Map(f = addDietCat,
               web = list("web1" = web1, "web2" = web1),
               missingDietcat = list("web1" = setdiff(dietcat, c(rowsWeb1, colsWeb1)),
                                     "web2" = setdiff(dietcat, c(rowsWeb2, colsWeb2))))
  web1 <- webLs$web1
  web1 <- webLs$web2

  ## calculate beta-diversity if networks are not empty
  temp <- suppressMessages(networkTempBetaDiv(web1 = web1,
                                              web2 = web2,
                                              dietcat = dietcat,
                                              method = method,
                                              mode = mode))
  temp[, PAGENAME := i]
  return(temp)
}


#' TAXONOMIC BETA-DIVERSITY FUNCTION
#' calculates beta diversity as the spp turnover between two networks.

#' @param web1 first adjency matrix (species network)
#' @param web2 second adjency matrix (species network)
#' @template dietcat
#' @param method is the method used to calculate beta diversity in 'composition' \code{mode}.
#'    Options are Sorensen's ("sorensen") and Simpson's ("simpson") beta-diversity
#'    and their nested component ("nested") (see Baselga 2010), the Jaccard
#'    dissimilarity index ("jaccard"), a multiplicative decomposition of alpha/gamma
#'    diversity (where alpha and gamma are calculated as inverse-Simpson; "decompSimpson"),
#'    all options ("all"), or any combination of them. All these options
#'    calculate beta-diversity (i.e., turnover) in terms of species composition.
#'    Defaults to "all"
#' @param mode can be "decompHills", "pairwiseHills" or "composition". "decompHills" and "pairwiseHills"
#'    use the \code{econetworks} package to calculate beta-diversity based on
#'    species and link composition, either through alpha/gamma decomposition
#'    ("decompHills") or by calculating pairwise beta-diversity ("pairwiseHills").
#'    In both cases, diversity is calculated using a Hills numbers approach. "composition" calculates
#'    beta-diversity based on species composition only, using one of \code{method}.
#'
#' @importFrom crayon blue
#' @importFrom data.table data.table
#' @importFrom econetwork divPartition disPairwise
#'
#' @export
networkTempBetaDiv <- function(web1, web2, dietcat, method = "all", mode = "composition") {
  if (!mode %in% c("decompHills", "pairwiseHills", "composition")) {
    stop("mode must be one of 'decompHills', 'pairwiseHills' or 'composition'")
  }

  ## make a named vector of output metrics
  betaMetrics <- c(sorensen = "beta_sore", simpson = "beta_simp",
                   nested = "beta_nest", jaccard = "diss_jacc",
                   decompSimpson = "beta_decompSimpson",
                   decompHills = "beta_decompHills",
                   pairwiseHills = "beta_pairwiseHills")

  if (mode == "composition") {
    method <- if (length(method) == 1 & "all" %in% method) {
      names(betaMetrics)
    }
    if (!all(method %in% names(betaMetrics)))
      stop("'method' must be 'all', or one or a combination of 'sorensen', 'simpson',
           'nested', 'jaccard', 'decompSimpson'")
    betaMetrics <- betaMetrics[method]
  } else {
    betaMetrics <- betaMetrics[mode]
  }

  if (all(is.na(web1)) & all(is.na(web2))) {
    temp <- data.table(PAGENAME = i,
                       beta_sore = NA, beta_simp = NA,
                       beta_nest = NA, diss_jacc = NA,
                       beta_decompSimpson = NA,
                       beta_decompHills = NA,
                       beta_pairwiseHills = NA)
  } else {
    if (mode == "composition") {
      ## convert to normal matrix
      if (!is(web1, "matrix")) {
        web1 <- as.matrix(web1)
      }
      if (!is(web2, "matrix")) {
        web2 <- as.matrix(web2)
      }

      ## make data.tables of spp presences per site
      ## create empty DTs first, because one of them can be empty
      ## this can happen even in BL, if climate introduced species that were not there previously - when calculating beta-div per family.
      blsite <- if (nrow(web1)) {
        data.table(SPP = rownames(web1), Pres = rep(1, nrow(web1)))
      } else {
        data.table(SPP = character(), Pres = numeric())
      }
      blsite <- blsite[!SPP %in% dietcat]

      futsite <- if (nrow(web2)) {
        data.table(SPP = rownames(web2), Pres = rep(1, nrow(web2)))
      } else {
        data.table(SPP = character(), Pres = numeric())
      }
      futsite <- futsite[!SPP %in% dietcat]

      ## make a site by spp matrix
      siteXspp <- merge(blsite, futsite, by = "SPP", all = TRUE, suffixes = c(".bl", ".fut"))

      ## convert to data.frame
      siteXspp <- as.data.frame(siteXspp)
      rownames(siteXspp) <- siteXspp$SPP     # make spp rownames
      siteXspp$SPP <- NULL
      siteXspp <- as.matrix(siteXspp[grepl("[[:alpha:]]", rownames(siteXspp)),])   ## convert to matrix and maintain only "true" spp (rownames can be "1" if a web had no spp)

      siteXspp[is.na(siteXspp)] <- 0    ## change NAs to absences
      siteXspp = t(siteXspp)            ## and transpose

      ## convert to data.frame again
      siteXspp <- as.data.frame(siteXspp)

      ## Computing Sorensen's and Simpson's beta-div + their nested component (Baselga 2010) +
      ## Jaccard dissimilarity index +
      ## beta as multiplicative decomposition of alpha/gamma div (calculated as inverse-Simpson)
      temp <- data.table(beta_sore = beta.sor(siteXspp)[1],
                         beta_simp = beta.sim(siteXspp)[1],
                         beta_nest = beta.nes(siteXspp)[1],
                         diss_jacc = vegdist(as.matrix(siteXspp), method = "jaccard", binary = TRUE)[1],
                         beta_decompSimpson = BetaDisQ(as.matrix(siteXspp), q = 0)["Pres.fut", "Pres.bl"],
                         beta_decompHills = NA,
                         beta_pairwiseHills = NA)
    } else {
      if (!is.null(method)) {
        message(blue("'method' is ignored in mode 'pairwiseHills' and 'decompHills'"))
      }

      ## if one network is empty, turnover is 1.
      if (all(is.na(web1)) | all(is.na(web2))) {
        temp <- data.table(PAGENAME = i,
                           beta_sore = NA, beta_simp = NA,
                           beta_nest = NA, diss_jacc = NA,
                           beta_decompSimpson = NA,
                           beta_decompHills = 1,
                           beta_pairwiseHills = 1)
      } else {
        pixelXspp.ls <- .convert2igraph(list("web1" = web1, "web2" = web2))

        if (mode == "pairwiseHills") {
          ## suppress progress bar
          invisible(capture.output(temp <- disPairwise(gList = pixelXspp.ls, type = "L")))
          temp <- unique(temp)
          temp <- data.table(beta_sore = NA,
                             beta_simp = NA,
                             beta_nest = NA,
                             diss_jacc = NA,
                             beta_decompSimpson = NA,
                             beta_decompHills = NA,
                             beta_pairwiseHills = temp)
        }

        if (mode == "decompHills") {
          temp <- divPartition(gList = pixelXspp.ls, type = "L", framework = "Chao")
          temp <- data.table(beta_sore = NA,
                             beta_simp = NA,
                             beta_nest = NA,
                             diss_jacc = NA,
                             beta_decompSimpson = NA,
                             beta_decompHills = temp$Beta,
                             beta_pairwiseHills = NA)
        }
      }
    }
  }
  temp <- temp[, ..betaMetrics]
  return(temp)
}


#' TEMPORAL BETA-DIVERSITY WRAPPER FUNCTION FOR PIXEL BY SPP MATRICES
#'
#' @description this function loads and prepares networks to
#'     calculate temporal spp turnover (beta diversity) and saves
#'     the outputs. Beta diversity is calculated per pairs of pixels
#'     of the same location, but different scenarios (a baseline
#'     scenario and a change scenario). This function was designed to be
#'     parallelized. The user should provide one scenario file at
#'     a time and a vector of baseline files (which will be matched to the
#'     scenario - SDM stat. model match)
#' @param method passed to \code{networkTempBetaDiv.master}
#' @param masterScen.file one scenario pixelXspp matrix file path, which
#'     can be made by \code{makeAndSaveMasterScen}.
#' @param masterBL.files a vectorlis of baseline pixelXspp matrix files, which
#'     can be made by \code{makeAndSaveMasterBL}.
#' @param scen.dir is the "root" directory where "scenario networks" were stored
#' @param out.dir is were  beta-diversity tables will be saved (in a folder tree
#'    respecting scenario folders.)
#' @param families character. determines for which families beta div will be calculated
#'    (this subsets the species). Can be all, one or a combination of
#'    \code{c("allTaxon", "bird", "mammal", "reptile", "amphibian")}, "allTaxon" meaning
#'    no subsetting is done and beta-diversity metrics are calculated across all
#'    species. Defaults to all of these.
#' @param toDo controls if beta-diversity hhas to be recalculated for 'all'
#'     scenarios or just 'missing' ones
#'
#' @return the output beta-diversity file path
#'
#' @importFrom tools file_path_sans_ext
#'
#' @export
calcTempBetaDiv.master <- function(masterScen.file, masterBL.files,
                                   families = c("allTaxon", "bird", "mammal", "reptile", "amphibian"),
                                   out.dir, scen.dir, toDo = "all", method = "all") {
  ## checks
  if (any(!families %in% c("allTaxon", "bird", "mammal", "reptile", "amphibian")))
    stop("'families' must be one, or several of
         c('allTaxon', 'bird', 'mammal', 'reptile', 'amphibian')")

  ## change "allTaxon" to nothing for file suffixes and subset
  ## selected taxon
  familiesSufs <- c("", "bird", "mammal", "reptile", "amphibian")
  names(familiesSufs) <- c("allTaxon", "bird", "mammal", "reptile", "amphibian")
  familiesSufs <- familiesSufs[families]

  cat(paste0("Start\n", date(), "\n"), append = TRUE)
  masterBL.file <- if (grepl("GAM", masterScen.file))
    grep("GAM", masterBL.files, value = TRUE) else
      grep("RF", masterBL.files, value = TRUE)

  ## make output file suffix and directory
  scen <- sub(scen.dir, "", dirname(masterScen.file))
  fileSuf <- file_path_sans_ext(basename(masterScen.file))
  fileSuf <- sub("masterdietFUND_", "", fileSuf)

  out.folder <- file.path(out.dir, scen)
  if (!dir.exists(out.folder))
    dir.create(out.folder, recursive = TRUE)

  ## make all file names and check for their existence
  ## only load webs if some files need to be created, or all computations need to be repeated
  outFiles <- file.path(out.folder, paste0("TaxonTempBetaDiv_", familiesSufs, "_", fileSuf, ".rds"))
  outFiles <- sub("__", "_", outFiles) ## remove double "__" when allTaxon is present
  names(outFiles) <- families
  rm(familiesSufs) ## clean-up to avoid confusion

  ## check which files need to be done, skip files done by subsetting families.
  if (any(file.exists(outFiles)) & toDo == "missing") {
    cat(paste0("Already done:\n",
               paste(outFiles[file.exists(outFiles)], collapse = "\n"),
               "\n... skipping\n"), append = TRUE)
    outFiles <- outFiles[!file.exists(outFiles)]
    families <- names(outFiles)
  }

  if (length(families)) {
    for (fam in families) {
      outFile <- outFiles[fam]
      charMatch <- if (fam != "allTaxon")
        paste0("^",  toupper(substr(fam, 1, 1)), "[[:xdigit:]]{1}") else NULL

      ## (re-)load masters
      masterBL <- tryCatch(readRDS(masterBL.file), error = function(e) NULL)
      masterScen <- tryCatch(readRDS(masterScen.file), error = function(e) NULL)

      if (exists("masterBL") & exists("masterScen")) {
        ## remove unnecessary columns
        rmCols <- grep("^(A|B|M|R)[[:digit:]]|PAGENAME", names(masterScen),
                       value = TRUE, invert = TRUE)
        colsBL <- setdiff(names(masterBL), rmCols)
        colsScen <- setdiff(names(masterScen), rmCols)

        ## subset taxon and "PAGENAME"
        if (!is.null(charMatch)) {
          colsBL <- grep(paste0("PAGENAME|(", charMatch, ")"), colsBL, value = TRUE)
          colsScen <- grep(paste0("PAGENAME|(", charMatch, ")"), colsScen, value = TRUE)
        }
        masterBL <- masterBL[, ..colsBL]
        masterScen <- masterScen[, ..colsScen]
        rm(colsBL, colsScen, rmCols)
        tempTaxonBetaDiv <- networkTempBetaDiv.master(pixXsppMat1 = masterBL,
                                                      pixXsppMat2 = masterScen,
                                                      method = "all")

        ## save data.table
        cat(paste0("saving...\n", outFile, "\n"), append = TRUE)
        saveRDS(tempTaxonBetaDiv, file = outFile)
        cat(paste0("End\n", date(), "\n"), append = TRUE)
      } else if (!exists("masterBL")) {
        cat(paste0("Loading failed for:\n", masterBL.file, "\n"), append = TRUE)
      } else if (!exists("masterScen")) {
        cat(paste0("Loading failed for:\n", masterScen.file, "\n"), append = TRUE)
      }
    }
  }
  return(outFiles)
}

#' TAXONOMIC BETA-DIVERSITY FUNCTION FOR PIXEL BY SPP MATRICES
#'
#' @description calculates beta diversity as the spp turnover between pairs of networks
#'    defined in two pixel by species matrices (containing species presences absences from
#'    networks in each pixel). Comparisons are made per pixel, between the two matrices,
#'    belonging from two different points in time, or two scenarios.
#'
#' @param pixXsppMat1 is a data.table of pixel id's (column "PAGENAME") and species (other columns)
#'    filled with 0s (absences from pixel network) and 1s (presences in pixel network)
#' @param pixXsppMat2 same as \code{pixXsppMat1}. The network compositions to be compared with
#'    \code{pixXsppMat1}
#' @param method is the method used to calculate beta diversity. Options are Sorensen's ("sorensen")
#'    and Simpson's ("simpson") beta-diversity and their nested component ("nested") (see Baselga 2010),
#'    the Jaccard dissimilarity index ("jaccard"), a multiplicative decomposition of alpha/gamma
#'    diversity (where alpha and gamma are calculated as inverse-Simpson; "decompSimpson"),
#'    all options ("all"), or any combination of them. Options "sorensen", "simpson", "nested", "jaccard",
#'    "decompSimpson" calculate species beta-diversity (i.e., turnover). Defaults to "all"
#'
#' @return a data.table of beta-diversity indices per pixel ID (column "PAGENAME")
#' @importFrom data.table rbindlist data.table melt
#'
#' @export
networkTempBetaDiv.master <- function(pixXsppMat1, pixXsppMat2,
                                      method = "all", ...) {
  ## make a named vector of output metrics
  betaMetrics <- c(sorensen = "beta_sore", simpson = "beta_simp",
                   nested = "beta_nest", jaccard = "diss_jacc",
                   decompSimpson = "beta_decompSimpson")
  method <- if (length(method) == 1 & "all" %in% method) {
    names(betaMetrics)
  }
  if (!all(method %in% names(betaMetrics))) {
    stop("'method' must be 'all', or one or a combination of 'sorensen', 'simpson',
           'nested', 'jaccard', 'decompSimpson'")
  }
  betaMetrics <- betaMetrics[method]

  ## save all unique pixel IDs for later
  pixIDs <- data.table(PAGENAME = union(pixXsppMat1$PAGENAME, pixXsppMat2$PAGENAME))

  ## remove pixels without species (they will be added back
  ## if they have species in one of the scenarios)
  cols <- setdiff(names(pixXsppMat1), "PAGENAME")
  pixXsppMat1 <- pixXsppMat1[pixXsppMat1[, rowSums(.SD), .SDcols = cols] > 0]
  cols <- setdiff(names(pixXsppMat2), "PAGENAME")
  pixXsppMat2 <- pixXsppMat2[pixXsppMat2[, rowSums(.SD), .SDcols = cols] > 0]

  pixXsppMat1[, Scen := "BL"]
  pixXsppMat2[, Scen := "Scen"]

  ## make sure each pixel exists for both scenarios
  pixIDScenCombs <- merge.data.table(pixXsppMat1[, .(PAGENAME, Scen)], pixXsppMat2[, .(PAGENAME, Scen)],
                                     by = "PAGENAME", all = TRUE)
  pixIDScenCombs[, `:=`(Scen.x = "BL",
                        Scen.y = "Scen")]  ## remove NAs
  pixIDScenCombs <- melt.data.table(pixIDScenCombs, id.vars = "PAGENAME", value.name = "Scen")

  pixXsppMat12 <- rbind(pixXsppMat1, pixXsppMat2, use.names = TRUE, fill = TRUE)

  ## join missing pixel/scenario combinations
  pixXsppMat12 <- pixXsppMat12[pixIDScenCombs[, .(PAGENAME, Scen)],
                               on = c("PAGENAME", "Scen"), nomatch = NA]
  if (any(is.na(pixXsppMat12)))
    pixXsppMat12 <- replaceNAs.data.table(pixXsppMat12, value = 0)

  rm(pixXsppMat1, pixXsppMat2)
  for (i in 1:10) gc()

  ## Computing Sorensen's and Simpson's beta-div + their nested component (Baselga 2010) +
  ## Jaccard dissimilarity index +
  ## beta as multiplicative decomposition of alpha/gamma div (calculated as inverse-Simpson)
  cols <- setdiff(names(pixXsppMat12), c("PAGENAME", "Scen"))
  setkey(pixXsppMat12, PAGENAME)

  tempTaxonBetaDiv <- pixXsppMat12[, .(beta_sore = beta.sor(.SD)[1],
                                       beta_simp = beta.sim(.SD)[1],
                                       beta_nest = beta.nes(.SD)[1],
                                       diss_jacc = vegdist(as.matrix(.SD), method = "jaccard", binary = TRUE)[1],
                                       beta_decompSimpson = BetaDisQ(as.matrix(.SD), q = 0)[2, 1]),
                                   .SDcols = cols, by = PAGENAME]
  cols <- c("PAGENAME", betaMetrics)
  tempTaxonBetaDiv <- tempTaxonBetaDiv[, ..cols]

  ## add back missing PAGENAMES with NAs
  tempTaxonBetaDiv[pixIDs, on = "PAGENAME", nomatch = NA]
  return(tempTaxonBetaDiv)
}



#' SPATIAL BETA-DIVERSITY WRAPPER FUNCTION
#'
#' This function loads and prepares networks to calculate temporal
#' spp turnover (beta diversity) and saves the outputs.
#' beta diversity is calculated per pairs of pixels of the same location,
#' but different scenarios (a baseline/current scenario and a change scenario).
#' this function was designed to be parallelized.
#'
#' @param ff is the a file path to list of networks obtained from a scenario of change.
#' @param mode can be either "decomp", to calculate beta-diversity as an
#'   alpha/gamma decomposition (using \code{econetwork::divPartition}), or "pairwise",
#'   to calculate pairwise beta-diversity between networks (using
#'   \code{econetwork::disPairwise}).
#' @param res.dir is the "root" directory where "scenario networks" were stored and out.dir is
#'   were results will be placed (in a folder tree respecting scenario folders.)
#' @param sampleNetworks controls if beta-diversity is to be calculated for a sample of networks (provide integer)
#'   or all networks (FALSE)
#' @param noReps determines number of reps for the sampling and calculating beta-diversity
#' @param networkGroups groups of networks within which to calculate beta-diversity. If provided,
#'   it should be a named list of vectors of network IDs. Sampling will occur within each group, so beware
#'   of sample size.
#' @param parallel activates parallel calculation of betadiv across samples. Uses future with plan(multiprocess)
#' @param noCores only used in 'parallel' is TRUE. defaults to 2.
#' @param toDo controls if beta-diversity has to be recalculated for 'all'
#'     scenarios or just 'missing' ones
#' @param cacheRepo passed to \code{reproducible::Cache}.
#'
#' @importFrom tools file_path_sans_ext
#' @importFrom reproducible Cache
#' @importFrom econetwork disPairwise divPartition
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#'
#' @export
calcSpaceBetaDiv <- function(ff, mode = "decomp", res.dir, out.dir, sampleNetworks = FALSE,
                             noReps = NULL, networkGroups = NULL, parallel = FALSE,
                             noCores = 2, toDo = "missing", cacheRepo = options("reproducible.cachePath")) {
  cat(paste0("Start\n", date(), ": ", ff, "\n"), append = TRUE)
  ## checks:
  if (!isFALSE(sampleNetworks) & is.null(noReps)) {
    warning("You chose to sample networks without repetition. We advise providing a 'noReps'")
    noReps <- 1
  }

  if (isFALSE(sampleNetworks) & !is.null(noReps)) {
    warning("You chose not to sample networks. Ignoring 'noReps'")
  }

  if (!isFALSE(sampleNetworks) & !is.numeric(sampleNetworks)) {
    stop("Provide a numeric/integer or FALSE to 'sampleNetworks'")
  }

  if (isFALSE(is.null(names(networkGroups)))) {
    if (!is.list(networkGroups)) {
      stop("'networkGroups' must be a named list")
    }
  }

  if (is.numeric(sampleNetworks) & isFALSE(is.null(networkGroups))) {
    if (any(sapply(networkGroups, function(x) length(x) < sampleNetworks))) {
      stop("Some groups in 'networkGroups' are smaller than 'sampleNetworks'")
    }
  }

  if (!mode %in% c("pairwise", "decomp")) {
    stop("mode must be one of 'pairwise' or 'decomp'")
  }

  ## make output file name and directory
  scen <- sub(res.dir, "", dirname(ff))
  fileSuf <- file_path_sans_ext(basename(ff))
  fileSuf <- sub("spp10kWdietFUND_", "", fileSuf)

  out.folder <- file.path(out.dir, scen)
  if (!dir.exists(out.folder))
    dir.create(out.folder, recursive = TRUE)

  outFile <- file.path(out.folder, paste0("SppLinkSpatialBetaDiv_", fileSuf, ".rds"))

  ## check if the file needs to be done, and skip if not
  if (file.exists(outFile) & toDo == "missing") {
    cat(paste0("Already done... skipping\n"), append = TRUE)
  } else {
    if (grepl("\\.R(D|d)ata", extension(basename(ff)))) {
      tryCatch(load(ff), error = function(e) NULL)
    } else if (grepl("\\.rds", extension(basename(ff)))) {
      tryCatch(pixelXspp.ls <- readRDS(ff), error = function(e) NULL)
    } else stop("ff must be a .RData/.Rdata/.rds file")

    if (!exists("pixelXspp.ls")) {
      cat(paste0("Loading failed for:\n", ff, "\n"), append = TRUE)
    } else {
      ## convert webs to igraph
      ## empty networks get NA
      ## note: cache the following lapply.
      pixelXspp.ls <- Cache(.convert2igraph,
                            networkList = pixelXspp.ls,
                            .cacheExtra = c(ff, length(pixelXspp.ls), object.size(pixelXspp.ls)),
                            userTags = c("pixelXspp2igraph"),
                            cacheRepo = cacheRepo,
                            omitArgs = c("userTags", "networkList"))
      # pixelXspp.ls <- loadFromCache(cacheId = "28bd7e698bf8e520") ##  NetworkSims/Baseline_SDMpres_GlobCover/ExtTrsh_quant10/spp10kWdietFUND_current_wm_bin_RF_noLUC_noIUCN_allhab_tresh0.rds
      ## remove missing webs
      NAwebs <- vapply(pixelXspp.ls, class, FUN.VALUE = character(1))
      pixelXspp.ls <- pixelXspp.ls[which(NAwebs == "igraph")]
      rm(NAwebs); for (i in 1:10) gc()

      if (isFALSE(is.null(networkGroups))) {
        noGroupNets <- setdiff(names(pixelXspp.ls), unlist(networkGroups))
        if (length(noGroupNets)) {
          warning(paste(length(noGroupNets), "networks were not found in 'networkGroups'",
                        "and will be ignored"))
        }
        noGroupNets <- setdiff(unlist(networkGroups), names(pixelXspp.ls))
        if (length(noGroupNets)) {
          warning(paste(length(noGroupNets), "pixels in 'networkGroups'",
                        "were not found in the pixel X network list and will be ignored"))
          networkGroups <- lapply(networkGroups, function(x, nets) x[x %in% nets],
                                  nets = names(pixelXspp.ls))
        }
      } else {
        networkGroups <- list("all" = names(pixelXspp.ls))
      }

      ## sample networks
      if (sampleNetworks) {
        message(paste("Sampling", sampleNetworks, "networks,", noReps, "times"))

        ## create sample names combining groups and reps (to avoid a nested list)
        sampNames <- unlist(lapply(names(networkGroups), function(x) paste0(x, "_samp", 1:noReps)))
        sampList <- lapply(sampNames, FUN = function(xx, x, size) {
          x <- x[[grep(sub("_.*", "", xx), names(x), value = TRUE)]]
          sample(x = x, size = size, replace = FALSE)
        }, x = networkGroups, size = sampleNetworks)
        names(sampList) <- sampNames
      } else {
        sampList <- networkGroups
        names(sampList) <- paste0(names(networkGroups), "_samp1")
      }

      if (mode == "pairwise") {
        ## SPATIAL TAXONOMIC BETA_DIVERSITY CALCULATION ---------------------------
        if (parallel) {
          plan(multisession, gc = TRUE, workers = 4)
          pairwiseBetaDivList <- future_lapply(sampList, FUN = function(samp, pixelXspp.ls, cacheRepo, outFile) {
            Cache(disPairwise,
                  gList = pixelXspp.ls[samp],
                  type = "L",
                  .cacheExtra = list(samp, object.size(pixelXspp.ls[samp]), outFile),
                  cacheRepo = cacheRepo,
                  omitArgs = c("gList"))
          },
          pixelXspp.ls = pixelXspp.ls,
          cacheRepo = cacheRepo,
          outFile = outFile)
          future:::ClusterRegistry("stop")
        } else {
          pairwiseBetaDivList <- lapply(sampList, FUN = function(samp, pixelXspp.ls, cacheRepo, outFile) {
            Cache(disPairwise,
                  gList = pixelXspp.ls[samp],
                  type = "L",
                  cacheRepo = cacheRepo,
                  .cacheExtra = list(samp, object.size(pixelXspp.ls[samp]), outFile),
                  omitArgs = c("gList"))
          },
          cacheRepo = cacheRepo,
          pixelXspp.ls = pixelXspp.ls,
          outFile = outFile)
        }
      } else {
        if (parallel) {
          plan(multisession, gc = TRUE, workers = 4)
          pairwiseBetaDivList <- future_lapply(sampList, FUN = function(samp, pixelXspp.ls, cacheRepo, outFile) {
            Cache(divPartition,
                  gList = pixelXspp.ls[samp],
                  type = "L",
                  framework = "Chao",
                  .cacheExtra = list(samp, object.size(pixelXspp.ls[samp]), outFile),
                  cacheRepo = cacheRepo,
                  omitArgs = c("gList"))
          },
          pixelXspp.ls = pixelXspp.ls,
          cacheRepo = cacheRepo,
          outFile = outFile)
          future:::ClusterRegistry("stop")
        } else {
          pairwiseBetaDivList <- lapply(sampList, FUN = function(samp, pixelXspp.ls, cacheRepo, outFile) {
            Cache(divPartition,
                  gList = pixelXspp.ls[samp],
                  type = "L",
                  framework = "Chao",
                  cacheRepo = cacheRepo,
                  .cacheExtra = list(samp, object.size(pixelXspp.ls[samp]), outFile),
                  omitArgs = c("gList"))
          },
          pixelXspp.ls = pixelXspp.ls,
          cacheRepo = cacheRepo,
          outFile = outFile)
        }
      }

      ## save data.table
      cat("saving...\n", append = TRUE)
      saveRDS(pairwiseBetaDivList, file = outFile)
      cat(paste0("End\n", date(), "\n"), append = TRUE)
    }
  }
}
