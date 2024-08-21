## ----------------------------------------------
## NETWORKS UTILITY FUNCTIONS
## ----------------------------------------------

#' CONVERT RASTER TO MATRIX OF PRESENCE/ABSENCES
#'
#' @param ras a `raster` with species presence data.
#'
#' @importFrom data.table data.table dcast
#' @export
ras2matrix <- function(ras) {
  ## ras is a raster file
  ## remove lakes
  ras[!is.na(lakes[])] <- NA

  ## convert to pixXhab matrics
  temp <- data.table(PAGENAME = as.character(mask10kID$PageName), Habs = ras[!is.na(mask10k[])])
  temp <- dcast(temp, formula = PAGENAME ~ Habs)
  temp2 <- as.data.frame(apply(temp[, !colnames(temp) %in% c("PAGENAME", "NA")],
                               MARGIN = c(1,2),
                               FUN = function(x) if (!is.na(x)) return(1) else(return(NA))))
  if (nrow(temp) == nrow(temp2)) {
    temp2$PAGENAME <- as.character(temp$PAGENAME)
  } else(stop("Problem with nrows"))

  ## PAGENAME as 1st column
  temp2 <- temp2[, c(ncol(temp2), 1:(ncol(temp2)-1))]

  return(temp2)
}


#' MAKE TABLE OF BIOGEOGRAPHIC REGION BY PIXEL ID
#'
#' Downloads and processes biogeographic regions layer
#' to identify biogeographi region names per network pixel ID
#'
#' @param rasterToMatch must be provided and match pixel IDs and names
#'   in `maskID`
#' @param maskID table with at least the column 'PageName' (the "names"
#'   of non-NA pixels). NA-ed pixels in rasterToMatch should not appear
#'   in this table.
#' @param args named list of other arguments passed to `reproducible::prepInputs`
#'   to download and process a layer of biogeographic regions. The
#'   following defaults are used, if not provided:
#'   `url = "https://www.eea.europa.eu/data-and-maps/data/biogeographical-regions-europe-3/zipped-shapefile-format-vector-polygon/zipped-shapefile-format-vector-polygon/at_download/file"`
#'   `archive = "file.zip"`
#'   `targetFile = "BiogeoRegions2016.shp"`
#'   `rasterToMatch = mask10k`
#'   `fun = "raster::shapefile"`
#'   `overwrite = TRUE`
#'   `destinationPath = "Habitats/Bioregions"`
#'   `useCache = FALSE`
#'   `cacheRepo = options("reproducible.cachePath")`
#'
#' @importFrom dplyr %>%
#' @importFrom data.table data.table setkey
#' @importFrom sf st_as_sf st_drop_geometry
#' @importFrom reproducible prepInputs Cache
#' @importFrom raster getValues
#' @export
makeBioregDT <- function(rasterToMatch = NULL, maskID = NULL, args = NULL) {
  if (any(is.null(rasterToMatch), is.null(maskID))) {
    stop("Please provide 'rasterToMatch' and 'maskID'")
  }

  if (sum(!is.na(rasterToMatch[])) != nrow(maskID)) {
    stop("nrow of 'maskID' must match number of non-NA pixels in 'rasterToMatch'")
  }

  if (!is.null(args) & any(!is.list(args), is.null(names(args)))) {
    stop("'args' must be a named list")
  }

  args$rasterToMatch <- rasterToMatch

  if (is.null(args$url)) {
    args$url <- paste0("https://www.eea.europa.eu/data-and-maps/data/biogeographical-regions-europe-3/",
                       "zipped-shapefile-format-vector-polygon/zipped-shapefile-format-vector-polygon/at_download/file")
  }
  if (is.null(args$archive)) {
    args$archive <- "file.zip"
  }
  if (is.null(args$targetFile)) {
    args$targetFile <- "BiogeoRegions2016.shp"
  }
  if (is.null(args$overwrite)) {
    args$overwrite <- TRUE
  }
  if (is.null(args$fun)) {
    args$fun <- "raster::shapefile"
  }
  if (is.null(args$destinationPath)) {
    args$destinationPath <- "Habitats/Bioregions"
  }

  if (is.null(args$cacheRepo)) {
    args$cacheRepo <- options("reproducible.cachePath")
  }

  bioregShp <- do.call(prepInputs, args = args)
  bioregShp <- st_as_sf(bioregShp)
  bioregShp$PK_UID <- as.numeric(bioregShp$PK_UID)

  ## rasterize - for now use gdal_rasterize, until getCover option is not available in fasterize
  bioregRas <- Cache(rasterizeCover,
                     rasterToMatch = rasterToMatch,
                     shp = bioregShp,
                     cacheRepo = args$cacheRepo,
                     field = "PK_UID",
                     noDataVal = 9999)

  ## "tabulate"
  bioregDT <- data.table(PAGENAME = maskID$PageName,
                         PK_UID = as.integer(getValues(bioregRas)[!is.na(getValues(bioregRas))]))

  ## merge with code names
  bioregCorresp <- bioregShp %>%
    st_drop_geometry(.) %>%
    data.table(.)
  bioregCorresp$PK_UID <- as.integer(bioregCorresp$PK_UID)
  setkey(bioregCorresp, PK_UID)
  setkey(bioregDT, PK_UID)

  ## some coastal pixels got "9999" - no bioregion
  bioregDT <- bioregCorresp[, .(PK_UID, code)][bioregDT]

  return(bioregDT)
}


#' FIND SPECIES IUCN STATUS
#' @param species character vector of species scientific names. see `rredlist::rl_search`
#' @param region region to search in. see `rredlist::rl_regions`
#' @param token IUCN API token. See https://apiv3.iucnredlist.org/api/v3/docs
#'   and `rredlist` documentation. Defaults to the token in this URL:
#'   https://apiv3.iucnredlist.org/api/v3/species/region/europe/page/0?token=9bb4facb6d23f48efbf424bb05c0c1ef1cf6f468393bc745d42179ac4aca5fee
#'
#' @export
findIUCNStatus <- function(species, region = NULL, token = "9bb4facb6d23f48efbf424bb05c0c1ef1cf6f468393bc745d42179ac4aca5fee") {
  if (!requireNamespace("rredlist", quietly = TRUE)) {
    stop("'rredlist' is not installed. Please install using:",
         "\ninstall.packages('rredlist')")
  }

  status <- sapply(species, FUN = function(x, region) {
    ## try checking if spp exists by searching in whole database
    if (length(rredlist::rl_search(x, key = token)$result)) {
      status <- rredlist:: rl_search(x, region = region, key = token)$result[["category"]]

      if (length(status)) {
        status
      } else {
        "NoAssessment"
      }
    } else {
      "Spp_not_found"
    }
  }, region = region)
  as.character(status)
}

#' FIND SPECIES ACCEPTED NAMES (IUCN)
#' @inheritParams findIUCNStatus
#'
#' @export
findIUCNAcceptedName <- function(species, token = "9bb4facb6d23f48efbf424bb05c0c1ef1cf6f468393bc745d42179ac4aca5fee") {
  if (!requireNamespace("rredlist", quietly = TRUE)) {
    stop("'rredlist' is not installed. Please install using:",
         "\ninstall.packages('rredlist')")
  }

  accptdName <- sapply(species, FUN = function(x){
    unique(rredlist::rl_synonyms(x, key = token)$result[["accepted_name"]])
  })
  accptdName
}


#' ADD MISSING DIET CATEGORIES (BASL NODES) TO A NETWORK
#'
#' @param web a square and named `matrix` representing an adjacency network
#' @param missingDietcat  the diet categories (i.e. basal nodes) that will
#'   be added to `web`
#'
#' @export
addDietCat <- function(web, missingDietcat) {
  temp <- cbind(matrix(0, nrow = nrow(web), ncol = length(missingDietcat)), web)
  temp <- rbind(matrix(0, nrow = length(missingDietcat), ncol = ncol(temp)), temp)
  rownames(temp) <- c(missingDietcat, rownames(web))
  colnames(temp) <- c(missingDietcat, colnames(web))
  return(temp)
}


#' REPLACE NAs IN DATA.TABLE  --------------------------------
#' replaces NAs in a data.table with a value
#'
#' @param DT is a data.table
#' @param value is the value to replace the NAs with. defaults to 0L
#'
#' @importFrom data.table set
#' @export
replaceNAs.data.table <- function(DT, value = 0L) {
  for (col in seq_along(DT))
    set(DT, i = which(is.na(DT[[col]])), j = col, value = value)
  DT
}


#' MAKE MAPS OF ROBUSTNESS METRICS
#'
#' @param ras is a raster of a robustness metric (like mean robustness across scenarios)
#' @param mask is an optional mask to be plotted under the robustness metric (i.e. to show the full study area as dark grey surface)
#' @param maskCol is the colour used to fill mask values. Defaults to `"grey30"`
#' @param factorLayer is a boolean indicating whether values in `ras` should be treated `as.factor`
#'
#' @importFrom ggplot2 ggplot geom_tile scale_fill_brewer scale_fill_distiller theme coord_equal
#' @importFrom ggpubr theme_pubr
#' @importFrom raster raster
#' @importFrom methods as
#' @export
robustnessMaps <- function(ras, mask = NULL, maskCol = NULL,
                           factorLayer = FALSE) {
  ## checks
  if (!is(ras, "RasterLayer"))
    stop("ras needs to be of class RasterLayer")
  if (!is.null(mask) & !is(mask, "RasterLayer"))
    stop("mask needs to be of class RasterLayer")

  if (is.null(maskCol))
    maskCol <- "grey30"

  ## make mask DF
  if (!is.null(mask)) {
    mask10k2 <- mask10k
    mask10k2[!is.na(mask10k2)] <- 1
    mask10k2 <- as.data.frame(as(mask10k2, "SpatialPixelsDataFrame"))
    names(mask10k2) <- c("layer", "x", "y")
    mask10k2 <- mask10k2[mask10k2$layer == 1, ]
  }

  plotData <- as.data.frame(as(ras, "SpatialPixelsDataFrame"))
  names(plotData) <- c("layer", "x", "y")

  plot1 <- ggplot(data = plotData)

  if (!is.null(mask)) {
    plot1 <- plot1 +
      geom_tile(data = mask10k2, aes(x = x, y = y), fill = maskCol, show.legend = FALSE)
  }

  plot1 <- if (factorLayer) {
    plot1 +
      geom_tile(aes(x = x, y = y, fill = as.factor(layer))) +
      scale_fill_brewer(palette = "RdYlBu", na.value = "transparent", direction = 1)
  } else {
    plot1 +
      geom_tile(aes(x = x, y = y, fill = layer)) +
      scale_fill_distiller(palette = "RdYlBu", na.value = "transparent", direction = 1)
  }

  plot1 <- plot1 +
    theme_pubr(base_size = 12, legend = "right") +
    theme(legend.key.height = unit(x = 0.07, units = "npc")) +
    coord_equal()

  return(plot1)
}



#' CONVERT LIST OF ADJENCY MATRICES TO IGRAPH
#'
#' Wrapper around `igraph::graph_from_adjacency_matrix` that
#'  outputs `NA` for empty networks.
#'
#' @param networkList is a list of adjacency matrices that can be converted to
#'   `igraph` using `igraph::graph_from_adjacency_matrix`
#'
#' @importFrom igraph graph_from_adjacency_matrix
.convert2igraph <- function(networkList) {
  lapply(networkList, FUN = function(web) {
    if (nrow(web) == 1) {
      return(NA)
    } else {
      graph_from_adjacency_matrix(adjmatrix = t(web), mode = "directed")
    }
  })
}
