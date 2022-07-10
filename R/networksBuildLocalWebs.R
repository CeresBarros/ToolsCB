## ----------------------------------------------
## FUNCTIONS TO SPATIALIZE NETWORKS
## (build local webs)
## Ceres May 2016, revamped >2018
## ----------------------------------------------

#' Function to calculate local webs - Baseline
#'
#' @template metaweb
#' @template SPPCODE
#' @template SPP.HAB
#' @template PIX.HAB
#' @template dietcat
#' @param HELP print help
#'
#' @export
BL_localweb <- function(metaweb = NULL, SPPCODE = NULL,
                        SPP.HAB = NULL, PIX.HAB = NULL,
                        dietcat = NULL, HELP = TRUE) {
  # Function to subset the meta-web square matrix according
  # a subgroup of species into a sub web square matrix
  if (HELP) warning("Description of HAB arguments: if spp x hab and pixel x hab matrices are provided habitat overlap will be used to filter spp links.
                   \n To avoid this message insert the following argument: HELP=FALSE")
  if (is.null(metaweb)) stop("Must provide a meta-web")
  if (is.null(SPPCODE)) stop("Must provide species codes")
  if (!is.null(SPP.HAB) & is.null(PIX.HAB)) stop("Must provide both spp x habitat and pixel x habitat matrices")
  if (is.null(SPP.HAB) & !is.null(PIX.HAB)) stop("Must provide both spp x habitat and pixel x habitat matrices")

  # -------------------------------------------------------
  # CALCULATE ORIGINAL LOCAL WEB
  # -------------------------------------------------------

  # calculate pixel potential web
  subweb <- metaweb[SPPCODE, SPPCODE]

  # FILTER 1. SPP REMOVAL BASED ON HABITAT OVERLAP
  # remove forbidden links, based on habitat overlap
  if (!is.null(SPP.HAB)) {
    sppcodes <- SPPCODE[!SPPCODE %in% dietcat]
    habs <- names(PIX.HAB[, -1])[PIX.HAB[, -1] > 0 & !is.na(PIX.HAB[, -1])]
    sppHAB <- as.matrix(SPP.HAB[sppcodes, habs, drop = FALSE])    # subset habitats by spp and habitats present in pixel
    sppHAB[is.na(sppHAB)] <- 0                           # NAs are considered 0s
    sppHAB[sppHAB > 1] <- 1                              # secondary and optimal habitats treated equally

    # adding diet categories with presences in all habitats (makes their addition to the final web easier)
    sppHAB <- rbind(matrix(1, nrow = length(dietcat), ncol = ncol(sppHAB), dimnames = list(dietcat, colnames(sppHAB))),
                    sppHAB)

    # creating a matrix of co-occurrences based on habitats
    sppXspp <- sppHAB %*% t(sppHAB)
    sppXspp[sppXspp > 1] <- 1                            # transforming to binary

    # multiplying each cell of the co-occurence matrix by the subweb will remove spp interactions that cannot occur in the pixel
    subweb <- sppXspp * subweb
  }


  # FILTER 2. REMOVAL BASED ON RESOURCES
  # remove SPP with no resources from metaweb (unconnected nodes + pure cannibals) iteratively
  while(any(rowSums(subweb[!row.names(subweb) %in% dietcat,, drop = FALSE], na.rm = TRUE) == 0)) {
    ## get spp that have prey
    SPPCODE <- c(dietcat, names(which(rowSums(subweb, na.rm = TRUE) > 0)))

    # update subweb, since some species were removed from predator list and have to be removed as prey
    subweb <- subweb[SPPCODE, SPPCODE]

    # remove pure cannibals
    pu.cannib <- colnames(subweb)[rowSums(subweb) == diag(subweb)]       # species (and diet categories) with sum of rows equal to diag represent pure cannibals
    pu.cannib2 <- pu.cannib[!pu.cannib %in% dietcat]                     # excluding diets from cannibal list

    subweb <- subweb[!rownames(subweb) %in% pu.cannib2, !colnames(subweb) %in% pu.cannib2] # remove pure cannibals
  }

  # remove unconnected nodes
  subweb <- subweb[which(!(rowSums(subweb, na.rm = TRUE) == 0 & colSums(subweb, na.rm = TRUE) == 0)), which(!(rowSums(subweb, na.rm = TRUE) == 0 & colSums(subweb, na.rm = TRUE) == 0)), drop = FALSE]

  return(subweb)
}


#' Function to build local networks according to habitat loss/LC change scenarios
#'
#' @template metaweb
#' @template SPPCODE
#' @template SPP.HAB
#' @template PIX.HAB
#' @template dietcat
#' @template EXT.TRSH
#' @template HELP
#'
#' @export
HabLoss_localweb <- function(metaweb = NULL, SPPCODE = NULL,
                             SPP.HAB = NULL, PIX.HAB = NULL,
                             EXT.TRSH = NULL, dietcat = NULL,
                             HELP = TRUE) {
  # Function to subset the meta-web square matrix according
  # a subgroup of species into a sub web square matrix
  if (HELP) warning("Description of HAB arguments: if spp x hab and pixel x hab matrices are provided habitat overlap will be used to filter spp links.
                   \n To avoid this message insert the following argument: HELP=FALSE")
  if (is.null(metaweb)) stop("Must provide a meta-web")
  if (is.null(SPPCODE)) stop("Must provide codes of species present in the BL web")
  if (!is.null(EXT.TRSH)) warning("Remove links between diet categories and non-basal predators")
  if (!is.null(SPP.HAB) & is.null(PIX.HAB)) stop("Must provide both spp x habitat and pixel x habitat matrices")
  if (is.null(SPP.HAB) & !is.null(PIX.HAB)) stop("Must provide both spp x habitat and pixel x habitat matrices")

  # -------------------------------------------------------
  # CALCULATE BASELINE LOCAL WEB
  # -------------------------------------------------------

  # calculate pixel BL web
  subweb <- metaweb[SPPCODE, SPPCODE]

  ## Add "prey" for diet categories
  init.prey <- rbind(as.matrix(EXT.TRSH), as.matrix(rowSums(subweb, na.rm = TRUE)[dietcat]))

  # -------------------------------------------------------
  # CALCULATE LOCAL WEB AFTER HABITAT CHANGES
  # -------------------------------------------------------

  # FILTER 1. SPP REMOVAL BASED ON HABITAT OVERLAP
  # PRIMARY EXTINCTIONS - REMOVAL OF FORBIDDEN LINKS
  # remove forbidden links, based on habitat overlap
  sppcodes <- SPPCODE[!SPPCODE %in% dietcat]
  habs <- names(PIX.HAB[, -1])[PIX.HAB[, -1] > 0 & !is.na(PIX.HAB[, -1])]
  sppHAB <- as.matrix(SPP.HAB[sppcodes, habs, drop = FALSE])    # subset habitats by spp and habitats present in pixel
  sppHAB[is.na(sppHAB)] <- 0                           # NAs are considered 0s
  sppHAB[sppHAB > 1] <- 1                              # secondary and optimal habitats treated equally

  # adding diet categories with presences in all habitats (makes their addition to the final web easier)
  sppHAB <- rbind(matrix(1, nrow = length(dietcat), ncol = ncol(sppHAB), dimnames = list(dietcat, colnames(sppHAB))),
                  sppHAB)

  # creating a matrix of co-occurrences based on habitats
  sppXspp <- sppHAB %*% t(sppHAB)
  sppXspp[sppXspp > 1] <- 1                            # transforming to binary

  # multiplying each cell of the co-occurence matrix by the subweb will remove spp interactions that cannot occur in the pixel
  subweb <- sppXspp * subweb

  ## get spp that after the habitat filter have some prey
  SPPCODE2 <- c(rownames(subweb)[rownames(subweb)%in% dietcat],
                rownames(subweb[rowSums(subweb, na.rm = TRUE) > 0,]))

  # FILTER 2. REMOVAL BASED ON RESOURCES
  # SECONDARY EXTINCTIONS
  # remove SPP with no resources from metaweb based on resource threshold (unconnected nodes + pure cannibals) iteratively
  while(any(rowSums(subweb, na.rm = TRUE) < init.prey[rownames(subweb),1])) {
    ## get spp that have enought prey
    SPPCODE3 <- names(which(rowSums(subweb, na.rm = TRUE) >= init.prey[rownames(subweb),1]))

    # update subweb, since some species were removed from predator list and have to be removed as prey
    subweb <- subweb[SPPCODE3, SPPCODE3]

    # remove pure cannibals
    pu.cannib <- colnames(subweb)[rowSums(subweb) == diag(subweb)]       # species (and diet categories) with sum of rows equal to diag represent pure cannibals
    pu.cannib2 <- pu.cannib[!pu.cannib %in% dietcat]                     # excluding diets from cannibal list

    subweb <- subweb[!rownames(subweb) %in% pu.cannib2, !colnames(subweb) %in% pu.cannib2] # remove pure cannibals
  }

  # remove unconnected nodes
  subweb <- subweb[which(!(rowSums(subweb, na.rm = TRUE) == 0 & colSums(subweb, na.rm = TRUE) == 0)), which(!(rowSums(subweb, na.rm = TRUE) == 0 & colSums(subweb, na.rm = TRUE) == 0)), drop = FALSE]

  # get present species after secondary extinctions
  SPPCODE3 <- row.names(subweb)

  ## -------------------------------------------------------
  ## CALCULATE NO. OF SPP EXTINCTIONS / INVASIONS
  ## -------------------------------------------------------

  ## remove diet categories first, for practical reasons
  SPPCODE <- SPPCODE[!SPPCODE %in% dietcat]
  SPPCODE2 <- SPPCODE2[!SPPCODE2 %in% dietcat]
  SPPCODE3 <- SPPCODE3[!SPPCODE3 %in% dietcat]

  ## Lists of primarily and secondarily extinct spp due to habitat changes
  Pext <- SPPCODE[!SPPCODE %in% SPPCODE2]
  Sext <- SPPCODE2[!SPPCODE2 %in% SPPCODE3]

  my_list <- list(web = subweb, Pext = Pext, Sext = Sext)
}

#' Function to build local networks according to spp removal scenarios
#'
#' @template metaweb
#' @template SPPCODE
#' @template SPP.HAB
#' @template PIX.HAB
#' @template dietcat
#' @param SPP.EXT vector of spp to be primarily extinct
#' @template EXT.TRSH
#' @template HELP
#'
#' @export
SppRm_localweb <- function(metaweb = NULL, SPPCODE = NULL,
                           SPP.HAB = NULL, PIX.HAB = NULL,
                           SPP.EXT = NULL, EXT.TRSH = NULL,
                           dietcat = NULL, HELP = TRUE) {
  # Function to subset the meta-web square matrix according
  # a subgroup of species into a sub web square matrix
  if (HELP) warning("Description of HAB arguments: if spp x hab and pixel x hab matrices are provided habitat overlap will be used to filter spp links.
                   \n To avoid this message insert the following argument: HELP=FALSE")
  if (is.null(metaweb)) stop("Must provide a meta-web")
  if (is.null(SPPCODE)) stop("Must provide codes of species present in the BL web")
  if (!is.null(SPP.HAB) & is.null(PIX.HAB)) stop("Must provide both spp x habitat and pixel x habitat matrices")
  if (is.null(SPP.HAB) & !is.null(PIX.HAB)) stop("Must provide both spp x habitat and pixel x habitat matrices")
  if (is.null(SPP.EXT)) stop("Must provide codes of primarily extinct species")

  # -------------------------------------------------------
  # CALCULATE LOCAL WEB AFTER PRIMARY EXTINCTIONS
  #                SECONDARY EXTINCTIONS
  # -------------------------------------------------------

  # remove primarily extinct species
  SPPCODE2 <- SPPCODE[!SPPCODE %in% SPP.EXT]

  # recalculate pixel web
  subweb <- metaweb[SPPCODE2, SPPCODE2]

  ## Add "prey" for diet categories to extinciton thresholds
  init.prey <- rbind(as.matrix(EXT.TRSH), as.matrix(rowSums(subweb, na.rm = TRUE)[dietcat]))

  # FILTER 1. REMOVAL BASED ON HABITAT OVERLAP
  # PRIMARY EXTINCTIONS - REMOVAL OF FORBIDDEN LINKS
  # remove forbidden links, based on habitat overlap
  if (!is.null(SPP.HAB)) {
    sppcodes <- SPPCODE[!SPPCODE %in% dietcat]
    habs <- names(PIX.HAB[, -1])[PIX.HAB[, -1] > 0 & !is.na(PIX.HAB[, -1])]
    sppHAB <- as.matrix(SPP.HAB[sppcodes, habs, drop = FALSE])    # subset habitats by spp and habitats present in pixel
    sppHAB[is.na(sppHAB)] <- 0                           # NAs are considered 0s
    sppHAB[sppHAB > 1] <- 1                              # secondary and optimal habitats treated equally

    # adding diet categories with presences in all habitats (makes their addition to the final web easier)
    sppHAB <- rbind(matrix(1, nrow = length(dietcat), ncol = ncol(sppHAB), dimnames = list(dietcat, colnames(sppHAB))),
                    sppHAB)

    # creating a matrix of co-occurrences based on habitats
    sppXspp <- sppHAB %*% t(sppHAB)
    sppXspp[sppXspp > 1] <- 1                            # transforming to binary

    # multiplying each cell of the co-occurence matrix by the subweb will remove spp interactions that cannot occur in the pixel
    subweb <- sppXspp * subweb
  }

  # FILTER 2. REMOVAL BASED ON RESOURCES
  # remove SPP with no resources from metaweb based on resource threshold (unconnected nodes) iteratively

  while(any(rowSums(subweb, na.rm = TRUE) < init.prey[rownames(subweb),1])) {
    ## get spp that have enought prey
    SPPCODE3 <- names(which(rowSums(subweb, na.rm = TRUE) >= init.prey[rownames(subweb),1]))

    # update subweb, since some species were removed from predator list and have to be removed as prey
    subweb <- subweb[SPPCODE3, SPPCODE3]

    # remove pure cannibals
    pu.cannib <- colnames(subweb)[rowSums(subweb) == diag(subweb)]       # species (and diet categories) with sum of rows equal to diag represent pure cannibals
    pu.cannib2 <- pu.cannib[!pu.cannib %in% dietcat]                       # excluding diets from cannibal list

    subweb <- subweb[!rownames(subweb) %in% pu.cannib2, !colnames(subweb) %in% pu.cannib2] # remove pure cannibals

  }

  # remove unconnected nodes
  subweb <- subweb[which(!(rowSums(subweb) == 0 & colSums(subweb) == 0)),which(!(rowSums(subweb) == 0 & colSums(subweb) == 0)), drop = FALSE]

  # update present species after secondary extinctions
  SPPCODE3 <- row.names(subweb)

  ## Lists of primarily and secondarily extinct spp due to habitat changes
  Pext <- SPPCODE[SPPCODE %in% SPP.EXT]
  Sext <- setdiff(SPPCODE2[!SPPCODE2 %in% dietcat], SPPCODE3[!SPPCODE3 %in% dietcat])

  my_list <- list(web = subweb, Pext = Pext, Sext = Sext)

  return(my_list)
}


#' Function to build local networks according to LC changes and changes in spp dists
#'
#' @template metaweb
#' @template SPP.HAB
#' @template PIX.HAB
#' @template ORIGSPP
#' @param PIX.HAB.ORIG baseline pixel x habitat `matrix`
#' @param PIX.HAB.FUT future pixel x habitat `matrix`
#' @param FUTSPP species IDs of future species presences
#' @template dietcat
#' @template EXT.TRSH
#' @template HELP
#'
#' @export
HabLoss_SppDist_localweb <- function(metaweb = NULL, ORIGSPP = NULL,
                                     SPP.HAB = NULL,
                                     PIX.HAB.ORIG = NULL, PIX.HAB.FUT = NULL,
                                     FUTSPP = NULL, EXT.TRSH = NULL,
                                     dietcat = NULL,
                                     HELP=TRUE) {
  ## Function to subset the meta-web square matrix according
  ## a subgroup of species into a sub web square matrix
  if (HELP) warning("Description of HAB arguments: if spp x hab and pixel x hab matrices are provided habitat overlap will be used to filter spp links.
                   \n To avoid this message insert the following argument: HELP=FALSE")
  if (is.null(metaweb)) stop("Must provide a meta-web")
  if (is.null(ORIGSPP) | is.null(FUTSPP)) stop("Must provide codes of species present in the baseline web and those of spp projected to be in the pixel in the future")
  if (!is.null(EXT.TRSH)) warning("Remove links between diet categories and non-basal predators")
  if ((!is.null(SPP.HAB) & (is.null(PIX.HAB.FUT) | is.null(PIX.HAB.ORIG))) |
      (is.null(SPP.HAB) & (!is.null(PIX.HAB.FUT) | !is.null(PIX.HAB.ORIG)))) {
    stop("Must provide both spp x habitat and future/current pixel x habitat matrices")
  }

  ## -------------------------------------------------------
  ## CALCULATE BASELINE LOCAL WEB
  ## -------------------------------------------------------

  ## calculate pixel web based on future spp distributions only
  subweb <- metaweb[FUTSPP, FUTSPP]

  ## Add "prey" for diet categories
  init.prey <- rbind(as.matrix(EXT.TRSH), as.matrix(rowSums(subweb, na.rm = TRUE)[dietcat]))

  ## -------------------------------------------------------
  ## CALCULATE LOCAL WEB AFTER HABITAT CHANGES
  ## -------------------------------------------------------

  ## FILTER 1. SPP REMOVAL BASED ON HABITAT OVERLAP
  ## SECONDARY EXTINCTIONS - REMOVAL OF FORBIDDEN LINKS
  ## remove forbidden links, based on habitat overlap
  ## here counted as secondary extinctions from spp distrib changes and/or habitat changes (excluding invasives)
  sppcodes <- FUTSPP[!FUTSPP %in% dietcat]
  habs <- names(PIX.HAB[, -1])[PIX.HAB[, -1] > 0 & !is.na(PIX.HAB[, -1])]
  sppHAB <- as.matrix(SPP.HAB[sppcodes, habs, drop = FALSE])    # subset habitats by spp and habitats present in pixel
  sppHAB[is.na(sppHAB)] <- 0                           # NAs are considered 0s
  sppHAB[sppHAB > 1] <- 1                              # secondary and optimal habitats treated equally

  ## adding diet categories with presences in all habitats (makes their addition to the final web easier)
  sppHAB <- rbind(matrix(1, nrow = length(dietcat), ncol = ncol(sppHAB), dimnames = list(dietcat, colnames(sppHAB))),
                  sppHAB)

  ## update spp that are present in future habitats (necessary in case false presences were not removed for spp dists)
  FUTSPP <- rownames(sppHAB)[rowSums(sppHAB) != 0]

  ## creating a matrix of co-occurrences based on habitats
  sppXspp <- sppHAB %*% t(sppHAB)
  sppXspp[sppXspp > 1] <- 1                            # transforming to binary

  ## multiplying each cell of the co-occurence matrix by the subweb will remove spp interactions that cannot occur in the pixel
  subweb <- sppXspp * subweb

  ## get spp that after the habitat filter have some prey
  FUTSPP2 <- c(rownames(subweb)[rownames(subweb)%in% dietcat],
               rownames(subweb[rowSums(subweb, na.rm = TRUE) > 0,]))

  # FILTER 2. REMOVAL BASED ON RESOURCES
  # SECONDARY EXTINCTIONS
  # remove SPP with no resources from metaweb based on resource threshold (unconnected nodes + pure cannibals) iteratively
  while(any(rowSums(subweb, na.rm = TRUE) < init.prey[rownames(subweb),1])) {
    ## get spp that have enought prey
    FUTSPP3 <- names(which(rowSums(subweb, na.rm = TRUE) >= init.prey[rownames(subweb),1]))

    # update subweb, since some species were removed from predator list and have to be removed as prey
    subweb <- subweb[FUTSPP3, FUTSPP3]

    # remove pure cannibals
    pu.cannib <- colnames(subweb)[rowSums(subweb) == diag(subweb)]       # species (and diet categories) with sum of rows equal to diag represent pure cannibals
    pu.cannib2 <- pu.cannib[!pu.cannib %in% dietcat]                     # excluding diets from cannibal list

    subweb <- subweb[!rownames(subweb) %in% pu.cannib2, !colnames(subweb) %in% pu.cannib2] # remove pure cannibals

  }

  # remove unconnected nodes
  subweb <- subweb[which(!(rowSums(subweb, na.rm = TRUE) == 0 & colSums(subweb, na.rm = TRUE) == 0)), which(!(rowSums(subweb, na.rm = TRUE) == 0 & colSums(subweb, na.rm = TRUE) == 0)), drop = FALSE]

  # get present species after secondary extinctions
  FUTSPP3 <- row.names(subweb)

  ## -------------------------------------------------------
  ## CALCULATE NO. OF SPP EXTINCTIONS / INVASIONS
  ## -------------------------------------------------------

  ## remove diet categories first, for practical reasons
  FUTSPP <- FUTSPP[!FUTSPP %in% dietcat]
  FUTSPP2 <- FUTSPP2[!FUTSPP2 %in% dietcat]
  FUTSPP3 <- FUTSPP3[!FUTSPP3 %in% dietcat]

  ## INVASIONS ---------------------------------------------
  ## initial invasions, needed for subsets - will updated later
  Invs <- FUTSPP[!FUTSPP %in% ORIGSPP]

  ## PRIMARY EXTINCTIONS -----------------------------------
  ## Primarily extinct spp (due to spp distrib changes and/or habitat changes)
  Pext <- ORIGSPP[!ORIGSPP %in% FUTSPP]

  ## SECONDARY EXTINCTIONS -----------------------------------

  ## Complex version
  # first round of secondary extinctions where habitat vs climate need to be distinguished
  Sext1 <- FUTSPP[!FUTSPP %in% FUTSPP2]
  Sext1 <- Sext1[!Sext1 %in% Invs]   # invasives that may not have settled do not count

  if (length(Sext1) > 0) {
    # Check what caused the primary extinction
    # build a web of ext + common species using original habitats and another using future habitats
    commonspp <- ORIGSPP[ORIGSPP %in% FUTSPP2]
    web.orighab <- BL_localweb(metaweb = metaweb, SPPCODE = c(dietcat, Sext1, commonspp),
                               SPP.HAB = SPP.HAB, PIX.HAB = PIX.HAB.ORIG,
                               dietcat = dietcat, HELP = FALSE)

    web.futhab <- BL_localweb(metaweb = metaweb, SPPCODE = c(dietcat, Sext1, commonspp),
                              SPP.HAB = SPP.HAB, PIX.HAB = PIX.HAB.FUT,
                              dietcat = dietcat, HELP = FALSE)

    # If spp did not have interactions with orig habitats, nor with fut habitats (so are excluded from the webs), then this was a climate secondary extinction (spp lost the resource due to climate -> remember that only COMMON spp and those with forbidden links removed are entering)
    Sext_sppdist <- intersect(Sext1[!Sext1 %in% rownames(web.orighab)], Sext1[!Sext1 %in% rownames(web.futhab)])

    # If spp whould have interactions with orig habitats (hence are present in web) but not in future (not present in fut web), then this was a habitat secondary extinction
    Sext_hab <- intersect(Sext1[Sext1 %in% rownames(web.orighab)], Sext1[!Sext1 %in% rownames(web.futhab)])

  } else {
    Sext_sppdist <- NULL
    Sext_hab <- NULL
  }

  # Now, had secondary extinctions coming from loss of resources after habitat filtering
  Sext2 <- FUTSPP2[!FUTSPP2 %in% FUTSPP3]
  Sext2 <- Sext2[!Sext2 %in% Invs]

  Sext_hab <- union(Sext_hab, Sext2)

  ## updated invasions
  Invs <- FUTSPP3[!FUTSPP3 %in% ORIGSPP]

  my_list <- list(web = subweb, Pext = Pext, Sext_sppdist = Sext_sppdist, Sext_hab = Sext_hab, Invs = Invs)
}


#' Function to build local networks according to changes in spp dists and, if desired, in habitats
#'
#' @template metaweb
#' @template SPP.HAB
#' @template PIX.HAB
#' @template ORIGSPP
#' @param FUTSPP species IDs of future species presences
#' @template dietcat
#' @template EXT.TRSH
#' @template HELP
#'
#' @export
SppDist_localweb <- function(metaweb = NULL, ORIGSPP = NULL,
                             SPP.HAB = NULL, PIX.HAB = NULL,
                             FUTSPP = NULL, EXT.TRSH = NULL,
                             dietcat = NULL,
                             HELP = TRUE) {
  ## Function to subset the meta-web square matrix according
  ## a subgroup of species into a sub web square matrix
  if (HELP) warning("Description of HAB arguments: if spp x hab and pixel x hab matrices are provided habitat overlap will be used to filter spp links.
                   \n To avoid this message insert the following argument: HELP = FALSE")
  if (is.null(metaweb)) stop("Must provide a meta-web")
  if (is.null(ORIGSPP) | is.null(FUTSPP)) stop("Must provide codes of species present in the baseline web and those of spp projected to be in the pixel in the future")
  if (!is.null(EXT.TRSH)) warning("Remove links between diet categories and non-basal predators")
  if (!is.null(SPP.HAB) & is.null(PIX.HAB)) stop("Must provide both spp x habitat and pixel x habitat matrices")
  if (is.null(SPP.HAB) & !is.null(PIX.HAB)) stop("Must provide both spp x habitat and pixel x habitat matrices")

  ## -------------------------------------------------------
  ## CALCULATE BASELINE LOCAL WEB
  ## -------------------------------------------------------

  ## calculate pixel web based on future spp distributions only
  subweb <- metaweb[FUTSPP, FUTSPP]

  ## Add "prey" for diet categories
  init.prey <- rbind(as.matrix(EXT.TRSH), as.matrix(rowSums(subweb, na.rm = TRUE)[dietcat]))

  ## -------------------------------------------------------
  ## CALCULATE LOCAL WEB AFTER HABITAT CHANGES
  ## -------------------------------------------------------

  ## FILTER 1. SPP REMOVAL BASED ON HABITAT OVERLAP
  ## SECONDARY EXTINCTIONS - REMOVAL OF FORBIDDEN LINKS
  ## remove forbidden links, based on habitat overlap
  ## here counted as secondary extinctions from spp distrib changes and/or habitat changes (excluding invasives)
  sppcodes <- FUTSPP[!FUTSPP %in% dietcat]
  habs <- names(PIX.HAB[, -1])[PIX.HAB[, -1] > 0 & !is.na(PIX.HAB[, -1])]
  sppHAB <- as.matrix(SPP.HAB[sppcodes, habs, drop = FALSE])    # subset habitats by spp and habitats present in pixel
  sppHAB[is.na(sppHAB)] <- 0                           # NAs are considered 0s
  sppHAB[sppHAB > 1] <- 1                              # secondary and optimal habitats treated equally

  ## adding diet categories with presences in all habitats (makes their addition to the final web easier)
  sppHAB <- rbind(matrix(1, nrow = length(dietcat), ncol = ncol(sppHAB), dimnames = list(dietcat, colnames(sppHAB))),
                  sppHAB)

  ## update spp that are present in future habitats (necessary in case false presences were not removed for spp dists)
  FUTSPP <- rownames(sppHAB)[rowSums(sppHAB) != 0]

  ## creating a matrix of co-occurrences based on habitats
  sppXspp <- sppHAB %*% t(sppHAB)
  sppXspp[sppXspp > 1] <- 1                            # transforming to binary

  ## multiplying each cell of the co-occurence matrix by the subweb will remove spp interactions that cannot occur in the pixel
  subweb <- sppXspp * subweb

  ## get spp that after the habitat filter have some prey
  FUTSPP2 <- c(rownames(subweb)[rownames(subweb)%in% dietcat],
               rownames(subweb[rowSums(subweb, na.rm = TRUE) > 0,]))

  ## FILTER 2. REMOVAL BASED ON RESOURCES
  ## SECONDARY EXTINCTIONS
  ## remove SPP with no resources from metaweb based on resource threshold (unconnected nodes + pure cannibals) iteratively
  while(any(rowSums(subweb, na.rm = TRUE) < init.prey[rownames(subweb),1])) {
    ## get spp that have enought prey
    FUTSPP3 <- names(which(rowSums(subweb, na.rm = TRUE) >= init.prey[rownames(subweb),1]))

    ## update codes and subweb, since some species were removed from predator list and have to be removed as prey
    subweb <- subweb[FUTSPP3, FUTSPP3]

    ## remove pure cannibals
    pu.cannib <- colnames(subweb)[rowSums(subweb) == diag(subweb)]       # species (and diet categories) with sum of rows equal to diag represent pure cannibals
    pu.cannib2 <- pu.cannib[!pu.cannib %in% dietcat]                     # excluding diets from cannibal list

    subweb <- subweb[!rownames(subweb) %in% pu.cannib2, !colnames(subweb) %in% pu.cannib2] # remove pure cannibals
  }

  ## remove unconnected nodes
  subweb <- subweb[which(!(rowSums(subweb, na.rm = TRUE) == 0 & colSums(subweb, na.rm = TRUE) == 0)), which(!(rowSums(subweb, na.rm = TRUE) == 0 & colSums(subweb, na.rm = TRUE) == 0)), drop = FALSE]

  ## get present species after secondary extinctions
  FUTSPP3 <- row.names(subweb)

  ## -------------------------------------------------------
  ## CALCULATE NO. OF SPP EXTINCTIONS / INVASIONS
  ## -------------------------------------------------------

  ## remove diet categories first, for practical reasons
  FUTSPP <- FUTSPP[!FUTSPP %in% dietcat]
  FUTSPP2 <- FUTSPP2[!FUTSPP2 %in% dietcat]
  FUTSPP3 <- FUTSPP3[!FUTSPP3 %in% dietcat]

  ## INVASIONS ---------------------------------------------
  ## initial invasions, needed for subsets - will updated later
  Invs <- FUTSPP[!FUTSPP %in% ORIGSPP]

  ## PRIMARY EXTINCTIONS -----------------------------------
  ## Primarily extinct spp (due to spp distrib changes and/or habitat changes)
  Pext <- ORIGSPP[!ORIGSPP %in% FUTSPP]

  ## SECONDARY EXTINCTIONS -----------------------------------
  ## Secondarily extinct spp
  Sext <- FUTSPP[!FUTSPP %in% FUTSPP3]
  Sext <- Sext[!Sext %in% Invs]   # invasives that may have settled do not count

  ## updated invasions
  Invs <- FUTSPP3[!FUTSPP3 %in% ORIGSPP]

  my_list <- list(web = subweb, Pext = Pext, Sext = Sext, Invs = Invs)
}
