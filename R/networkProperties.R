#' ----------------------------------------------
#' NETWORKS PROPERTIES/ANALYSIS FUNCTIONS
#' ----------------------------------------------

#' CALCULATE MEAN BETA DIVERSITY (PAIRWISE DISTANCE
#' Summarizes network pairwise distances into mean distance from a
#'   focal network to others
#'
#' @param outputs a data.table with columns \code{PAGENAME}  (the network ID)
#'   and \code{meanDistance}
#' @param distObj a \code{dist} object calculated using \code{econetwork::disPairwise}
#'
#' @importFrom data.table data.table
#'
#' @export
calcNetMeanDistances <- function(distObj) {
  ## convert to matrix
  tempMatrix <- as.matrix(distObj)
  ## remove diagonal to calculate network mean distance to others
  diag(tempMatrix) <- NA

  ## calculate mean distances and create data.table
  meanDistances <- rowMeans(tempMatrix, na.rm = TRUE)
  sdDistances <- apply(tempMatrix, 1, sd, na.rm = TRUE)
  meanDistances <- data.table(PAGENAME = names(meanDistances), meanDistance = meanDistances,
                              sdDistance = sdDistances)
  return(meanDistances)
}


#' CALCULATE MIN PREY EXTINCTION THRESHOLDS
#'
#' @param x is a list of file paths of to the list of networks (full)
#' @param quants defines the quantile values to use
#' @param out.dir is the destination directory for outputs
#' @template dietcat
#'
#' @importFrom plyr rbind.fill
#'
#' @export
calc.extThresh <- function(x, quants, out.dir, dietcat) {
  ## loading all BL networks
  if (grepl("\\.R(D|d)ata", extension(basename(x)))) {
    load(x)
  } else if (grepl("\\.rds", extension(basename(x)))) {
    pixelXspp.ls <- readRDS(x)
  } else stop("x must be a .RData/.Rdata/.rds file")
  outfile.suff <- paste0(sub("base.*", "base",
                             sub(".*current", "current",
                                 file_path_sans_ext(basename(x)))),
                         ".txt")
  ## getting the number of prey per spp in each pixel
  pixXprey <- rbind.fill(lapply(1:length(pixelXspp.ls), FUN = function(i){
    pix = names(pixelXspp.ls[i])
    web = as.matrix(pixelXspp.ls[[i]])
    web = web[!rownames(web) %in% dietcat, ] ## remove diet categories as predators

    if(sum(dim(web)) > 0){
      prey <- as.data.frame(t(as.matrix(rowSums(web))))
      rownames(prey) = pix
      return(prey)
    }
  }))


  ## calculating quantiles and median values per species
  ext.thresh <- do.call(rbind.data.frame, apply(as.matrix(pixXprey), 2, FUN = function(x, quants){
    temp <- data.frame(min = min(x, na.rm = TRUE), median = median(x, na.rm = TRUE))
    colnames(temp) = c("min", "median")

    temp2 <- data.frame(as.list(quantile(x, probs = quants, na.rm = TRUE)))
    colnames(temp2) = sub("\\.", "_", paste0("quant", quants*100))

    temp3 <- cbind(temp, temp2)

    return(temp3)
  }, quants = quants))

  ## save outputs
  write.table(pixXprey, file = file.path(out.dir, paste0("Number_preyXpixel_BLwebs_", outfile.suff)))
  write.table(ext.thresh, file = file.path(out.dir,paste0("Ext_thresh_BLwebs_", outfile.suff)))
}



#' CALCULATE NETWORK METRICS
#' Function to calculate several network metrics
#'
#' @param web a square \code{matrix} representing an adjacency matrix.
#' @param normalise use normalised version of generality, vulnerability (and their
#'   standard deviations)?
#'
#' @importFrom igraph graph.adjacency walktrap.community modularity degree
#'
#' @export
netw.metrics <- function(web = NULL, normalise = FALSE, verbose = TRUE) {
  if (is.null(web)) stop("Must provide a network matrix prey x predator")

  ## no. species, no. links and connectance
  S <- length(web[,1])
  L <- sum(web)
  C <- L/(S^2)

  ## Modularity
  igraph.mat <- graph.adjacency(web, mode = "directed")    # converting web to igraph object
  comm <- walktrap.community(igraph.mat)                 # calculating communities
  Q <- modularity(igraph.mat, membership = membership(comm))    # calculating modularity

  ## Parameters of degree distribution
  deg <- degree(graph.adjacency(web))            # calculate degree distribution (basically calculates the degree)
  # fit power law, can give an error if there are unconnected spp (basal consumers that aren't predated)
  degdist <- try(fitting2(deg), silent = FALSE)
  if (any(grepl("Error", degdist))) {
    dd.alpha <- NA
    dd.beta <- NA
  } else {
    dd.alpha <- degdist$alpha                      # get parameters of the truncated power law
    dd.beta <- degdist$beta
  }

  ## Plot observed and fitted degree distributions.
  # gamma.trace <- 1 - pgamma((0:length(Pk)), shape = degdist$alpha, scale = degdist$beta)
  # plot(log(1:(length(Pk))),log(rev(cumsum(rev(Pk)))),
  #      pch=3, xlab="log(k)", ylab="log(cumulative P(k))")
  # lines(log(1:(length(Pk)+1)), log(gamma.trace), lty=1, lwd=2, col = "red")

  ## Taxa - proportion of basal, intermediate and top consumers, and omnivorous spp
  propB <- sum(colSums(web)==0 )/ nrow(web)
  propT <- sum(rowSums(web)==0 & colSums(web)!=0 )/ nrow(web)
  propI <- 1- (propB + propT)

  if (normalise) {
    if (verbose) warning("Calculating normalised metrics of generality/vulnerability")
    normGen <- MeanGenerality_norm(web)
    normVul <- MeanVulnerability_norm(web)
    SDnormGen <- SDGenerality_norm2(web)
    SDnormVul <- SDVulnerability_norm2(web)

    ## join metrics
    tab <- as.matrix(data.frame(S = S, L = L, C = C, Q = Q,
                                dd.alpha = dd.alpha, dd.beta = dd.beta,
                                normGen = normGen, normVul = normVul, SDnormGen = SDnormGen, SDnormVul = SDnormVul,
                                propB = propB, propI = propI, propT = propT))
  } else {
    Gen <- MeanGenerality(web)
    Vul <- MeanVulnerability(web)
    SDGen <- SDGenerality(web)
    SDVul <-  SDVulnerability(web)

    ## join metrics
    tab <- as.matrix(data.frame(S = S, L = L, C = C, Q = Q,
                                dd.alpha = dd.alpha, dd.beta = dd.beta,
                                Gen = Gen, Vul = Vul, SDGen = SDGen, SDVul = SDVul,
                                propB = propB, propI = propI, propT = propT))
  }

  return(tab)
}


#' CALCULATE TROPHIC LEVEL STATISTICS
#' Function to calculate trophic level stats
#'
#' @param web a square \code{matrix} representing an adjacency matrix.
#' @template dietcat
#'
#' @export
tlstats <- function(web = NULL, dietcat = dietcat) {
  ## Use PredationMatrixToLinks() to create a Cheddar community from a predation
  community <- Community(data.frame(node = colnames(web)), trophic.links = PredationMatrixToLinks(web),
                         properties = list(title = "BL"))

  community <- RemoveCannibalisticLinks(community, title='community');

  ## get all spp trophic levels
  tl <- PreyAveragedTrophicLevel(community)

  ## Omnivory
  resource.spp <- ResourcesByNode(community) ## get resource spp for each predator
  n.resources <- sapply(resource.spp, length) ## get no. resources

  resource.tl <- sapply(resource.spp, FUN = function(x) {
    return(length(unique(tl[x])))   ## get the number of different trophic levels predated upon
  })

  omnivs <- n.resources >= 2 & resource.tl >= 2  ## a spp is an omnivore if it predates 2 or more spp of different trophic levels

  if (!is.null(dietcat)) {
    omnivs <- omnivs[!names(omnivs) %in% dietcat]   ## if there are DCs exclude them
  }

  tab <- as.matrix(data.frame(propOmn = sum(omnivs)/length(omnivs),   ## proportion of omnivores (omnivory)
                              mean.TL = mean(tl), ## mean trophic level
                              max.TL = max(tl),   ## max trophic level
                              sd.TL = sd(tl)))    ## sd trophic level
  return(tab)
}


#' DETECT INTERACTIONS FUNCTION
#'
#' Function to detect potential prey/predators
#'
#' @template metaweb
#' @template SPPCODE
#' @param MODE one of "SUBSET", "EATS" or "EATEN". Use \code{HELP = TRUE} for a
#'   description of each option
#' @param HELP print description of \code{MODE}?
#'
#' @export
potential.inter <- function(metaweb = NULL, SPPCODE = NULL, MODE = NULL, HELP = TRUE) {
  if (HELP) warning("Description of MODE:
                   \n    SUBSET selects a subset of the metaweb based on a group of species (SPPCODE);
                   \n    EATS selects all species that a species eats;
                   \n    EATEN selects all species that eat a species.
                   \n To avoid this message insert the following argument: HELP=FALSE")
  if (is.null(MODE)) stop("Provide MODE argument - SUBSET, EATS, EATEN")
  if (length(SPPCODE) > 1 && MODE=="EATS") stop("EATS only works for single species")
  if (length(SPPCODE) > 1 && MODE=="EATEN") stop("EATEN only works for single species")

  if (MODE=="EATS")   subweb <- metaweb[SPPCODE, metaweb[SPPCODE,, drop = FALSE] >=1, drop =FALSE]
  if (MODE=="EATEN")  subweb <- metaweb[metaweb[,SPPCODE, drop = FALSE] >= 1, SPPCODE, drop = FALSE]

  return(subweb)
}


#' SPECIES CENTRALITY FUNCTION
#'
#' Function to calculate the centrality of a list of spp within a network, or a list of networks
#' network can be a single network or a list of networks to extract spp centrality values from
#' if network is a list, then species also needs to be a list of equal length
#'
#' @param networkList a named \code{list} of networks, even if containing a single network
#' @param sppExclude list of species to exclude from degree centrality calculations
#'   (all other metrics are calcualted on the full network)
#' @param metric centrality metric to calculate. Either "degree"
#'   (using \code{igraph::degree(network, mode = "all", v = spp)}),
#'   "betweenness" (using \code{igraph::betweenness(network, directed = TRUE)})
#'   "eigenvector" centrality (using \code{eigen_centrality(graph = network, directed = FALSE, scale = FALSE)}).
#'   (see Bauer et al 2010, Ecological Complexity).
#' @param ... passed to reproducible::Cache
#'
#' @importFrom reproducible Cache
#'
#' @export
spp.centrality <- function(networkList, sppExclude, metric = c("degree", "betweenness", "eigenvector"),
                           ...) {
  ## do some data checks
  if (is.list(networkList)) {
    if (is.null(names(networkList))) {
      stop("network must have `names()`")
    }
  } else {
    stop("network must be a list")
  }

  args <- list(metric = metric)
  if (!missing(sppExclude)) {
    args$sppExclude <- sppExclude
  }

  centrality.ls <- Cache(Map,
                         f = .calcCentrality,
                         network = networkList,
                         MoreArgs = args,
                         ...)
  ## to test for errors.
  # centrality.ls <- Map(f = function(pix, networkList, sppExclude, metric) {
  #   print(pix)
  #   .calcCentrality(network = networkList[["AA469"]], sppExclude = sppExclude,
  #                   metric = metric)
  # }, pix = names(networkList),
  # MoreArgs = list(networkList = networkList,
  #                 sppExclude = sppExclude,
  #                 metric = metric))

  return(centrality.ls)
}


#' CALCULATE SPECIES CENTRALITY
#'
#' Internal function to calculate centrality of a list of species in a network
#'
#' @param network a network coercible to \code{matrix} and \code{igraph}
#' @param sppExclude list of species to exclude the network before calculating degree
#'   centrality (all other metrics are calculated on the full network)
#' @param metric centrality metric to calculate. Either "degree"
#'   (using \code{igraph::degree(network, mode = "all", v = spp)}),
#'   "betweenness" (using \code{igraph::betweenness(network, directed = TRUE)})
#'   "eigenvector" centrality (using \code{eigen_centrality(graph = network, directed = FALSE, scale = FALSE)}).
#'   (see Bauer et al 2010, Ecological Complexity).
#'
#' @importFrom igraph graph_from_adjacency_matrix degree eigen_centrality betweenness
#' @importFrom data.table data.table
.calcCentrality <- function(network, sppExclude,
                            metric = c("degree", "betweenness", "eigenvector")) {
  ## get network
  network <- as.matrix(network)

  if (isFALSE(all(is.na(network)))) {
    if (!identical(sort(rownames(network)), sort(colnames(network)))) {
      stop("network is not and adjency network matrix, column and row names do not match")
    }

    # sppNames <- grep("PAGENAME|x", names(master.ext), value = TRUE, invert = TRUE)
    # ## get spp list and extract spp that exist in pix
    # spp <- master.ext[PAGENAME == pix, ..sppNames]
    # spp <- names(spp)[spp == 1]
    spp <- colnames(network)

    ## ---------------------------------------------
    ## CALCULATING CENTRALITY MEASURES -------------
    ## I chose degree and betweenness centrality (see Bauer et al 2010, Ecological Complexity)
    if (length(spp)) {
      ## Converting adjancency to igraph format
      network <- graph_from_adjacency_matrix(adjmatrix = t(network))

      ## Degree
      if ("degree" %in% metric) {
        if (!missing(sppExclude)) {
          spp2 <- setdiff(spp, sppExclude)
        } else {
          spp2 <- spp
        }

        network2 <- graph_from_adjacency_matrix(adjmatrix = t(network[spp2, spp2, drop = FALSE]))
        degr <- degree(network2, mode = "all")
        if (any(degr == 0)) {
          warning(paste("There are unconnected nodes in this network"))
        }
        centrality <- data.table(degr = degr, Spp = names(degr))
      }

      ## Vertex/node betweenness
      if ("betweenness" %in% metric) {
        VB <- betweenness(graph = network, directed = TRUE)
        VB <- data.table(betw = VB, Spp = names(VB))

        if (exists("centrality")) {
          centrality <- centrality[VB, on = "Spp"]
        } else {
          centrality <- VB
        }
      }

      if ("eigenvector" %in% metric) {
        eigenC <- eigen_centrality(graph = network, directed = FALSE, scale = FALSE)
        eigenC <- data.table(eigen = eigenC$vector, Spp = names(eigenC$vector))
        if (exists("centrality")) {
          centrality <- centrality[eigenC, on = "Spp"]
        } else {
          centrality <- eigenC
        }
      }
    }
  }

  if (!exists("centrality", inherits = FALSE)) {
    centrality <- NA
  }

  return(centrality)
}


#' NETWORK METRICS PCA FUNCTION
#'
#' Calculates a PCA on a matrix of network metrics usin\code{ade4::dudi.pca}.
#'   The function automatically excludes covariates with >0.9 correlation
#'
#' @param metricsDT is a data.table of metrics (columns) calculated for each network (rows)
#' @param PLOT logical. controls whether a PCA biplot is produced.
#' @param dim is used for faster caching and should be dim(metricsDT)
#'
#' @importFrom ade4 dudi.pca
#' @import ggplot2
#'
#' @export
calcNetworkBLMetricsPCA <- function(metricsDT, PLOT = FALSE, dim) {
  cols <- setdiff(names(metricsDT), "PAGENAME")
  ## find covariates that are highly correlated
  corMectricsDT <- cor(metricsDT[, ..cols], use = "complete.obs")
  ## focus on lower triangulaFr part and exclude diagonal
  corMectricsDT[!lower.tri(corMectricsDT)] <- NA
  highCorInds <- which(corMectricsDT > 0.9, arr.ind = TRUE)
  ## remove "self-correlations"
  highCorInds <- highCorInds[which(rowSums(highCorInds) != highCorInds[,1] * 2),]

  ## select columsn to keep for PCA
  cols <- setdiff(rownames(corMectricsDT), row.names(highCorInds))

  ## PCA
  metricsDT <- na.omit(metricsDT)

  cols <- setdiff(names(metricsDT), c("PAGENAME"))
  metricsPCA <- dudi.pca(metricsDT[, ..cols], nf = length(cols), scannf = FALSE, center = TRUE)

  if (PLOT) {
    # Plot axes loadings
    plot1 <- ggplot() +
      geom_point(data = metricsPCA$li,
                 mapping = aes(x = Axis1, y = Axis2),
                 color = "black", size = 1) +
      geom_segment(data = metricsPCA$c1,
                   mapping = aes(x = 0, y = 0, xend = CS1, yend = CS2),
                   arrow = arrow(length = unit(0.2, "cm")), alpha = 0.75, color = "red") +
      geom_text(data = metricsPCA$c1,
                mapping = aes(x = CS1, y = CS2, label = cols), size = 5, vjust = 1.2, color = "red") +
      geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
      scale_x_continuous(limits = c(-1,1)) +
      scale_y_continuous(limits = c(-1,1)) +
      theme_bw() + labs(x = "PC 1", y = "PC 2")

    print(plot1)
  }

  metricsPCA$tab$PAGENAME <- metricsDT$PAGENAME
  return(metricsPCA)
}


#' AVERAGING NETWORK METRICS FUNCTION
#'
#' Calculates a set of summary statistical metrics on a data.table of
#' network metrics per pixel, per group of a grouping variable (or several)
#' @param DT is a data.table of variables on which the summary stats will be calculated
#' @param stats is a list of summary statistics to calculate. Select from "mean", "sd", "cv",
#'  "median", "min", "max", where "cv" is the coefficient of variation
#' @param id.vars character vector of grouping variable names
#' @param measure.vars character vector variables to summarise. If NULL then all non id.vars are used
#' @param na.rm is passed to summary statistics functions
#'
#' @export
calcNetworkMetricsSummStats <- function(DT, stats = c("mean", "sd", "cv", "median", "min", "max"),
                                        id.vars, measure.vars = NULL, na.rm = FALSE) {
  ## Checks
  if (!any(class(DT) == "data.table"))
    DT <- as.data.table(DT)
  if (is.null(id.vars)) {
    stop("must provide at least one variable in id.vars")
  } else {
    if (!any(names(DT) %in% id.vars))
      stop(paste("Can't find", setdiff(id.vars, names(DT)), "in DT"))
  }

  if (is.null(measure.vars)) {
    measure.vars <- setdiff(names(DT), id.vars)
  } else {
    if (!any(names(DT) %in% measure.vars))
      stop(paste("Can't find", setdiff(measure.vars, names(DT)), "in DT"))
  }
  if (any(!stats %in% c("mean", "sd", "cv", "median", "min", "max")))
    stop("stats must be one, or several, of 'mean', 'sd', 'cv', 'median', 'min', 'max'")

  if ("median" %in% stats)
    medians <- DT[, lapply(.SD, median, na.rm = na.rm), by = id.vars, .SDcols = measure.vars]
  if ("mean" %in% stats)
    means <- DT[, lapply(.SD, mean, na.rm = na.rm), by = id.vars, .SDcols = measure.vars]
  if ("sd" %in% stats)
    sds <- DT[, lapply(.SD, sd, na.rm = na.rm), by = id.vars, .SDcols = measure.vars]

  if (all(c("mean", "sd", "cv") %in% stats)) {
    cvs <- sds[, ..measure.vars] / means[, ..measure.vars]
    cvs <- cbind(sds[, ..id.vars], cvs)
  } else {
    ## it's faster to calculate mean/sd separately first
    means <- DT[, lapply(.SD, mean, na.rm = na.rm), by = id.vars, .SDcols = measure.vars]
    sds <- DT[, lapply(.SD, sd, na.rm = na.rm), by = id.vars, .SDcols = measure.vars]

    cvs <- sds[, ..measure.vars] / means[, ..measure.vars]
    cvs <- cbind(sds[, ..id.vars], cvs)
  }

  if ("min" %in% stats)
    mins <- DT[, lapply(.SD, min, na.rm = na.rm), by = id.vars, .SDcols = measure.vars]
  if ("max" %in% stats)
    maxs <- DT[, lapply(.SD, max, na.rm = na.rm), by = id.vars, .SDcols = measure.vars]

  ## merge data.tables
  objList <- lapply(grep(paste(paste0("^", stats, "s$"), collapse = "|"),
                         ls(), value = TRUE),
                    FUN = function(x) get(x))

  names(objList) <- grep(paste(paste0("^", stats, "s$"), collapse = "|"),
                         ls(), value = TRUE)

  ## change variable names
  objList <- lapply(names(objList), FUN = function(x) {
    DT <- objList[[x]]

    oldNames <- grep(paste(id.vars, collapse = "|"), names(DT), invert = TRUE, value = TRUE)
    setnames(DT, oldNames, paste0(oldNames, "_", sub("s$", "", x)))

    DT
  })

  ## join all tables
  objList <- lapply(objList, FUN = function (DT) setkeyv(DT, id.vars))
  my.merge <- function(x,y) merge(x,y)
  summDT <- Reduce(my.merge, objList)

  summDT
}
