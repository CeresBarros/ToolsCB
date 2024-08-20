#############################################################################################
#                                                                                           #
#     Alpha, beta, gamma decomposition with parametrization of the dominance (q) effect     #
#                                                                                           #
#     Refs : Chalmandrier et al. Ecology                                                    #
#            Chao et al. 2010                                                               #
#            Leinster & Cobbold 2012                                                        #
#                                                                                           #
#############################################################################################

#' Site diversity according to Leinster & Cobbold 2012.
#'
#' `divLeinster` calculates the diversity of each site of a site by species
#'   matrix according to the q parameter (which provide more or less weight to
#'   species abundance), following Leinster & Cobbold 2012.
#'
#' @param spxp sites (row) by species (cols) `matrix` with or without `rownames` and `colnames`
#' @param Z similarity `matrix` used into all functions.
#' @param q parameter influencing the weight of species abundances on
#'   diversity calcualtions as per Leinster & Cobbold 2012.
#' @param check  should arguments be checked?
#' @export
divLeinster <- function(spxp, Z = NULL, q = 2, check = TRUE){
  #Calcul the diversity of each site of sites by species matrix.
  #spxp columns and Z rows and columns are assumed to be in the same order.
  if (is.null(Z)) Z <- diag(ncol(spxp))
  if (check){
    if (!inherits(spxp, "matrix")) {
      stop("object \"spxp\" is not of class \"matrix\"")}
    if (!inherits(Z, "matrix")) {
      stop("object \"Z\" is not of class \"matrix\"")}
    if (!all(c(ncol(Z), nrow(Z)) == ncol(spxp))){
      stop("object \"Z\" and object \"spxp\" does not have matching dimensions")}
  }
  spxp <- sweep(spxp, 1, rowSums(spxp), "/")
  Zp <- Z %*% t(spxp)

  if (q != 1 & q != Inf){
    mat <- t(spxp) * (Zp)^(q-1)
    mat[is.na(mat)] <- 0
    D <- colSums(mat) ^ (1/(1-q))
  }
  if (q == Inf)  {
    D <- 1 / apply(Zp, 2, max)
  }
  if (q == 1){
    D <- apply(Zp^t(spxp), 2, function(x) 1/prod(x))
  }
  return(D)
}

#' alpha-, beta- and gamma-diversity multiplicative decomposition
#'
#' `abgDecompQ`  performs a alpha, beta, gamma multiplicative decomposition
#'   using Leinster's diversity indices.
#'
#' @inheritParams divLeinster
#' @param mult determines whether beta-diversity  should be calculated using
#'   multiplicative decomposition (the default) or as a percentage
#'   (B = (G-mA)/G * 100), where mA is mean alpha-diversity
#'
#' @export
abgDecompQ <- function(spxp, Z = NULL, q = 2, mult = TRUE, check = TRUE) {
  #Calcul the diversity of each site of sites by species matrix.
  #spxp columns and Z rows/cols are assumed to be in the same order.
  if (is.null(Z)) Z <- diag(ncol(spxp))
  if (check){
    if (!inherits(spxp, "matrix")) {
      stop("object \"spxp\" is not of class \"matrix\"")}
    if (!inherits(Z, "matrix")) {
      stop("object \"Z\" is not of class \"matrix\"")}
    if (!all(c(ncol(Z), nrow(Z)) == ncol(spxp))){
      stop("object \"Z\" and object \"spxp\" does not have matching dimensions")}
  }

  site.weight <- rep(1/nrow(spxp), nrow(spxp))
  spxp <- sweep(spxp, 1, rowSums(spxp), "/")

  gamma.ab <- colSums(sweep(spxp, 1, site.weight, "*"), na.rm = TRUE)

  Gamma <- divLeinster(t(as.matrix(gamma.ab)), Z=Z , q=q, check = FALSE)
  Alphas <- divLeinster(spxp, Z=Z , q=q, check = FALSE)

  if (q != 1 & q != Inf) {
    mAlpha <- sum(site.weight * (Alphas ^ (1 - q)), na.rm = TRUE)^(1 / (1 - q))
  }
  if (q==1){
    mAlpha <- exp(sum(site.weight * log(Alphas), na.rm = TRUE))
  }
  if (q==Inf){
    mAlpha <- min(Alphas)
  }
  if (mult==TRUE)
  {
    Beta <- Gamma / mAlpha
  } else(Beta <- ((Gamma - mAlpha) / Gamma) * 100)    # added by Ceres

  names(Alphas) <- row.names(spxp)
  res <- list(Gamma=Gamma, Beta=Beta, mAlpha=mAlpha, Alphas=Alphas)

  return(res)
}

#' Pairwise beta-diversity calculation
#'
#' `BetaDisQ` calculates the pairwise beta-diversity (minus 1) between sites of a site by species matrix according to the q parameter using the afformentionned functions
#' Allows a parametrization of the dominance effect
#'
#' @inheritParams divLeinster
#' @inheritParams abgDecompQ
#' @export
BetaDisQ <- function(spxp, Z = NULL, q = 2, check = TRUE, mult = TRUE){
  #Calcul the site pairwise diversity of a sites by species matrix.
  #spxp columns and Z rows/cols are assumed to be in the same order.
  if (is.null(Z)) Z <- diag(ncol(spxp))
  if (check){
    if (!inherits(spxp, "matrix")) {
      stop("object \"spxp\" is not of class \"matrix\"")}
    if (!inherits(Z, "matrix")) {
      stop("object \"Z\" is not of class \"matrix\"")}
    if (!all(c(ncol(Z), nrow(Z)) == ncol(spxp))){
      stop("object \"Z\" and object \"spxp\" does not have matching dimensions")}
  }

  N <- nrow(spxp)
  dis <- matrix(NA, N, N)
  for (i in 2:N) {
    for (j in 1:(i-1)) {
      spxp.dummy <- spxp[c(i,j), ]
      res <- abgDecompQ(as.matrix(spxp.dummy), Z = Z, q = q, check = check, mult = mult)
      dis[i, j] <- dis[j, i] <- res$Beta
    }
  }

  diag(dis) <- 1
  dis <- dis - 1
  row.names(dis) <- colnames(dis) <- row.names(spxp)
  return(dis)
}

#' Data prep for diversity metrics
#'
#' `chaoObjects` is a data preparation function.
#'   It returns adequate arguments for `abgDecompQ`, `BetaDisQ` and
#'   `divLeinster` to perform a diversity analysis using Chao's diversity index.
#'   Warning: this formula only works with an ultrametric tree!
#'
#' @param spxp sites (row) by species (cols) `matrix` with or without `rownames` and `colnames`
#' @param phy ultrametric phylogenetic tree.

chaoObjects <- function(spxp, phy) {

  if (!"ape" %in% installed.packages()) {
    stop("Please install ape package")
  } else {
    requireNamespace("ape")
  }

  if (!"phangorn" %in% installed.packages()) {
    stop("Please install phangorn package")
  } else {
    requireNamespace("phangorn")
  }

  if (!inherits(phy, "phylo")){
    stop("object \"phy\" is not of class \"phylo\"")}
  if (!inherits(spxp, "matrix")) {
    stop("object \"spxp\" is not of class \"matrix\"")}
  if (ncol(spxp) != length(phy$tip.label)){
    stop("object \"phy\" and object \"spxp\" does not have the same number of species")}

  Ancestors.sp <- lapply(1:length(phy$tip.label), function(x) c(Ancestors(phy, x, "all"), x))
  Branches <- lapply(Ancestors.sp, function(x) which((phy$edge[,1] %in% x) & (phy$edge[,2] %in% x)))
  Li <- unlist(lapply(Branches, function(x) sum(phy$edge.length[x]))) #Tip - root distances

  ultra <- all.equal.numeric(var(unlist(Li)), 0, tolerance = 1e-7) #Is it ultrametric
  if (!ultra) stop ("object \"phy\" must be an ultrametric tree")

  freq.dummy <- unlist(lapply(1:length(Branches),function(i) rep(i,length(Branches[[i]]))))
  desc <- lapply(unlist(Branches), function(x) Descendants(phy, phy$edge[x,2], type ="tips")[[1]])

  tmp <- lapply(desc, function(desc_i){
    x <- rep(0, length(freq.dummy))
    x[freq.dummy %in% desc_i] <- 1
    return(x)
  })
  Z <- do.call(rbind, tmp)
  pi <- sweep(spxp[,freq.dummy], 2, phy$edge.length[unlist(Branches)], FUN = "*")/Li[1]

  res <- list(Z = Z,
              pi = pi)
  return(res)
}
