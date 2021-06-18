#' CLEAN SF OBJECT FROM DUPLICATED COLUMNS
#'
#' this functions removes potentially duplicated data columns in an sf object
#'
#' @param sfObj is an sf object from the sf package
#'
#' @export
#'
#' @return an \code{sf} object without duplicated columns
#'
#' @importFrom sf st_set_geometry st_geometry st_geometry<-

sfRmDupNACols <- function(sfObj) {
  transposed <- t(st_set_geometry(sfObj, NULL))
  ## get duplicates and NAs (note that cols became rows)
  dupCols <- duplicated(transposed)
  NAcols <- sapply(sfObj[,, drop = TRUE], FUN = function(var) all(is.na(var)))

  if (any(dupCols)) {
    transposed <- transposed[!dupCols,]
    dataSF <- data.frame(t(transposed), stringsAsFactors = FALSE)

    ## convert columns to numeric where appropriate
    numCols <- intersect(names(dataSF), names(which(sapply(st_set_geometry(sfObj, NULL), is.numeric))))
    for (col in numCols) {
      dataSF[, col] <- as.numeric(dataSF[, col])
    }

    ## check for NA cols that were not removed
    NAcols <- sapply(dataSF, FUN = function(var) all(is.na(var)))
    if (any(NAcols))
      dataSF <- dataSF[, !NAcols]

    ## add geometry to data.frame to re-make sf object
    st_geometry(dataSF) <- st_geometry(sfObj)
    sfObj <- dataSF
  } else {
    if (any(NAcols)) {
      sfObj <- sfObj[, !NAcols]
    } else message("no duplicated or NA columns were found.")
  }
  return(sfObj)
}


#' RENAME AND SUBSET FIELDS IN SF OBJECT ACCORDING TO TABLE
#'
#' @inheritParams sfRmDupNACols
#' @param namesTable is a two-column \code{data.frame} with the names to be replaced (1st column) and the new names (2nd column)
#'  if there are NAs in the second columns, the original columns will be removed if \code{rmNAs == TRUE}
#' @param rmNAs determines whether columns with missing new names are removed or not. Defaults to \code{TRUE}
#'
#' @export
#'
#' @return an \code{sf} object
#'
#' @importFrom sf st_set_geometry st_geometry

renameCleanSfFields <- function(sfObj, namesTable, rmNAs = TRUE) {
  origNames <- names(st_set_geometry(sfObj, NULL))
  if (!all(origNames %in% namesTable$Name_shp))
    stop("Some names in 'sfObj' are missing from 'namesTable'")

  ## get new names
  rownames(namesTable) <- namesTable[, 1]
  newNames <- as.character(namesTable[origNames, 2])

  ## if there are missing new names, either remove them, or maintain original name
  if (any(is.na(newNames))) {
    if (rmNAs) {
      toRM <-  origNames[is.na(newNames)]
      sfObj[, toRM] <- NULL

      ## re-do steps above
      origNames <- names(st_set_geometry(sfObj, NULL))
      newNames <- as.character(namesTable[origNames, 2])
    } else {
      newNames[is.na(newNames)] <- origNames[is.na(newNames)]
    }
  }

  names(sfObj)[names(sfObj) %in% origNames] <- newNames

  return(sfObj)
}


#' VALIDATE GEOMETRIES
#'
#' Checks for corrupt geometries in an \code{sf} object
#' and validates them using \code{sf::st_make_valid}

#' @inheritParams sfRmDupNACols
#' @param dim is used for faster caching
#'
#' @export
#'
#' @return an \code{sf} object with validated geometries. Or
#'   or an error if geometries could not be validated.
#'
#' @importFrom sf st_is_valid st_make_valid
#' @importFrom stats na.omit

validateGeomsSf <- function(sfObj, dim) {
  ## Any corrupt or invalid geometries?
  if (any(is.na(st_is_valid(sfObj))) |
      any(na.omit(st_is_valid(sfObj)) == FALSE)) {
    message("Invalid gemoetries found. Attempting to make valid using st_make_valid")
    sfObj <- st_make_valid(sfObj)
  }

  if (any(is.na(st_is_valid(sfObj))) |
      any(na.omit(st_is_valid(sfObj)) == FALSE))
    stop("st_make_valid did not work, please check what's wrong.")

  return(sfObj)
}


#' DOWNLOAD KMZ AND CONVERT TO SHAPEFILE
#'
#' @param url is a Google Drive URl
#' @param archive is the .zip file name where the .kmz is contained
#' @param destinationPath path is the path to where archive will downloaded to, and where the final shapefile will be saved.
#' @param overwrite passed to \code{googledrive::drive_download} and \code{sp::shapefile}
#'
#' @export
#'
#' @return a \code{shapefile}
#'
#' @importFrom googledrive drive_download as_id
#' @importFrom sf st_read st_zm as_Spatial
#' @importFrom raster shapefile
#' @importFrom utils unzip

prepKMZ2shapefile <- function(url, archive, destinationPath, overwrite = TRUE) {
  ## check
  if (class(url) != "character" | is.null(url))
    stop("Provide url as a character string")
  if (!grepl("\\.zip$", archive))
    stop("archive should be the name of the .zip file in url")
  if (is.null(destinationPath)) {
    message("archive will be downloaded to tempdir()")
    destinationPath <- tempdir()
  }

  archivePath <- file.path(destinationPath, archive)
  downloaded_file <- drive_download(file = as_id(url),
                                    path = archivePath,
                                    overwrite = overwrite)
  if (downloaded_file$local_path != archivePath |
      downloaded_file$name != archive)
    stop(paste0("Downloaded file name (", downloaded_file$name, ") and archive (",
                archive, ") do not match."))

  ## unzip .zip and .kmz
  fileKMZ <- unzip(archivePath, exdir = destinationPath)
  fileKML <- unzip(fileKMZ, exdir = destinationPath)

  ## load as sf
  sfObj <- st_read(fileKML)
  if (any(names(sfObj) %in% "Description"))
    sfObj$Description <- NULL  ## weird unnecessary column
  ## convert to shapefile
  shpObj <- as_Spatial(st_zm(sfObj))

  ## delete unecessary .kmz/.kml files and save .shp
  file.remove(fileKML, fileKMZ)
  shpFile <- sub(".zip", ".shp", archive, fixed = TRUE)
  shapefile(shpObj, filename = file.path(destinationPath, shpFile), overwrite = overwrite)

  return(shpObj)
}

#' DRAW CONVEX HULL AROUND POLYGON
#'
#' Draws a convex hull around vertice points of a polygon \code{shapefile}.
#'
#' @param x a \code{SpatialPolygons}, or a \code{SpatialPolygonsDF} object
#' NOTE: this function as been passed to amc.
#'
#' @export
#'
#' @return a convex hull polygon
#'
#' @importFrom sp SpatialPoints polygons

outerBuffer <- function(x) {
  if (class(x) == "SpatialPolygons" | class(x) == "SpatialPolygonsDataFrame") {
    ## Get polygon vertices
    pts <- SpatialPoints(do.call(rbind, lapply(x@polygons, FUN = function(x) {
      return(x@Polygons[[1]]@coords)
    })))

    ## Draw convex hull around points and extract polygons slot
    hull <- polygons(convHull(pts))

    return(hull)
  } else(stop("x must be a SpatialPolygons, or SpatialPolygonsDF"))
}

#' JOIN SPATIAL OBJECTS FUNCTION
#'
#' joins \code{shapefiles} or \code{rasters}
#'
#' @param files is a character string of file names
#' @param destinationPath passed to \code{reproducible::prepInputs}
#' @param urls a character vector or list of URLs for each file
#'
#' @export
#'
#' @return joined \code{shapefile} or \code{raster}
#'
#' @importFrom reproducible prepInputs
#' @importFrom raster crs bind

loadBindSpatialObjs <- function(files, destinationPath, urls = NULL) {
  ## name URLs with file names
  if (length(urls) > 1) {
    names(urls) <- files

    if (all(grepl(".shp", files))) {
      spObj.ls <- lapply(files, FUN = function(targetFile) {
        prepInputs(targetFile = file.path(destinationPath, targetFile),
                   url = urls[targetFile], destinationPath = destinationPath,
                   fun = "shapefile", pkg = "raster")
      })
      if (length(unique(sapply(spObj.ls, FUN = function(spObj) as.character(crs(spObj))))) == 1) {
        joined = do.call(bind, spObj.ls)
      } else(stop("Files do not share the same projection"))

    }  else {
      ## check if all are raster files
      if (all(grepl(".grd", files)) |
         all(grepl(".asc", files)) |
         all(grepl(".tif", files)) |
         all(grepl(".img", files))) {

        spObj.ls <- lapply(files, FUN = function(targetFile) {
          prepInputs(targetFile = file.path(destinationPath, targetFile),
                     url = urls[targetFile], destinationPath = destinationPath,
                     fun = "raster", pkg = "raster")
        })

        ## check projections match and do the join
        if (length(unique(sapply(spObj.ls, FUN = function(spObj) as.character(crs(spObj))))) == 1) {
          joined = do.call(bind, spObj.ls)
        } else(stop("Files do not share the same projection"))

      } else stop("All files should be in the same format (.shp, .grd, .asc, .tif or .img)")
    }

  } else {
    if (all(grepl(".shp", files))) {
      spObj.ls <- lapply(files, FUN = function(targetFile) {
        prepInputs(targetFile = targetFile,
                   url = urls, destinationPath = destinationPath,
                   fun = "shapefile", pkg = "raster")
      })
      if (length(unique(sapply(spObj.ls, FUN = function(spObj) as.character(crs(spObj))))) == 1) {
        joined = do.call(bind, spObj.ls)
      } else(stop("Files do not share the same projection"))

    }  else {
      ## check if all are raster files
      if (all(grepl(".grd", files)) |
         all(grepl(".asc", files)) |
         all(grepl(".tif", files)) |
         all(grepl(".img", files))) {

        spObj.ls <- lapply(files, FUN = function(targetFile) {
          prepInputs(targetFile = targetFile,
                     url = urls, destinationPath = destinationPath,
                     fun = "raster", pkg = "raster")
        })

        ## check projections match and do the join
        if (length(unique(sapply(spObj.ls, FUN = function(spObj) as.character(crs(spObj))))) == 1) {
          joined = do.call(bind, spObj.ls)
        } else(stop("Files do not share the same projection"))

      } else stop("All files should be in the same format (.shp, .grd, .asc, .tif or .img)")
    }
  }
  return(joined)
}


#' FIND NEIGHBOURS IN MATRIX
#'
#' finds the 8 neighbours of each cell in a matrix
#'
#' @param mat a \code{matrix}. Border of NA's is added to the matrix
#'
#' @export
#'
#' @return a matrix of 8 rows, with a columns per cell of the input matrix,
#'  which is treated by columns

neighboursMatrix <- function(mat) {
  if (!is(mat, "matrix")) {
    stop("mat needs to be a 'matrix'")
  }
  mat2 <- cbind(NA, rbind(NA, mat ,NA), NA)   ## makes a border of NAs
  addresses <- expand.grid(x = 1:nrow(mat), y = 1:ncol(mat)) ## all matrix coordinates
  neighs <- c()
  for (i in 1:-1) {
    for (j in 1:-1) {
      if (i != 0 || j != 0) {
        neighs <- rbind(neighs,
                        mat2[addresses$x + i + 1 + nrow(mat2)*(addresses$y + j)])   ## each column contains the neighbours a cell (going by columns in mat)
      }
    }
  }
  return(neighs)
}

## CHECK PROJECTIONS ------------------------------------
## sfObj.list is a list of spatial objects
checkProjections <- function(sfObj.list){
  projs <- sapply(sfObj.list, FUN = function(x) {
    return(
      eval(expr = parse(text = paste0("projection(", x, ")")))
    )
  })
  return(projs)
}


#' CROP & MASK TO STUDY AREA
#'
#' @param study.area is a \code{Raster} or \code{Spatial} object
#' @param tocrop a \code{Raster}
#' @param method passed to \code{raster::projectRaster} - might need to be changed for factors
#'
#' @export
#'
#' @importFrom raster projectRaster crs crop mask

cropToStudyArea <- function(study.area, tocrop, method = "bilinear") {
  temp <- tocrop

  if (class(study.area) == class(temp) & class(study.area) == "RasterLayer") {
    temp <- projectRaster(from = temp, to = study.area, crs = crs(study.area), method = method)
  } else {
    if (crs(study.area)@projargs != crs(tocrop)@projargs){
      temp <- projectRaster(tocrop, crs = crs(study.area))
    }
  }
  temp <- crop(x = temp, y = study.area)
  temp <- mask(x = temp, mask = study.area)

  return(temp)
}

#' VECTOR TO BINARY MATRIX
#'
#' converts a vector of values into a binary "presence/absence" matrix
#'
#' @param x is a vector. Values must be coercible to character
#'
#' @export
#'
#' @return a presence/absence \code{matrix}
#'
#' @importFrom stats model.matrix
vector2binmatrix <- function(x) {
  x <- as.character(x)
  return(model.matrix( ~ x-1))
}