
################## BEGIN NODE METRICS ##################

#' In-degree per node
#' @template M
InDegree <- TrophicGenerality <- NumberOfResources <- function(M){
  return(colSums(M));
}

#' Out-degree per node
#' @template M
OutDegree <- TrophicVulnerability <- NumberOfCosumers <- function(M){
  return(rowSums(M));
}

#' Node degree
#' @template M
Degree <- function(M){
  return(InDegree(M)+OutDegree(M));
}

#' Node normalised Generality
#' @template M
NormalisedGenerality <- function(M){
  return(TrophicGenerality(M)/(sum(M)/dim(M)[1]));
}

#' Node normalised Vulnerability
#' @template M
NormalisedVulnerability <- function(M){
  return(TrophicVulnerability(M)/(sum(M)/dim(M)[1]));
}

#' Calculate if a spp is an omnivore based on the number of different trophic levels predated by a spp
#' @template M
#' @param level a function to calculate the trophic level. See \code{cheddar}
#'   for options
IsOmnivore2 <- function(M, level = PreyAveragedTrophicLevel){
  #Use PredationMatrixToLinks() to create a Cheddar community from a predation
  community <- Community(data.frame(node = colnames(M)), trophic.links = PredationMatrixToLinks(M),
                         properties = list(title = "Test2"))
  community <- RemoveCannibalisticLinks(community, title='community');

  # get resource spp for each predator
  resource.spp <- ResourcesByNode(community)

  #get no. resources
  n.resources <- sapply(resource.spp, length)

  # get all spp trophic levels
  tl <- level(community)

  # get the number of different trophic levels predated upon
  resource.tl <- sapply(resource.spp, FUN = function(x){
    return(length(unique(tl[x])))
  })

  # a spp is an omnivore if it predates 2 or more spp of different trophic levels
  return(n.resources >= 2 & resource.tl >= 2)
}

################## END NODE METRICS ##################

################## BEGIN NETWORK METRICS ##################

#' Network mean Generality
#' @template M
MeanGenerality <- function(M){
  return(sum(colSums(M))/sum((colSums(M)!=0)));
}

#' Network mean Vulnerability
#' @template M
MeanVulnerability <- function(M){
  return(sum(rowSums(M))/sum((rowSums(M)!=0)));
}

#' Network SD Generality
#' @template M
SDGenerality <- function(M){
  return(sd(InDegree(M)[InDegree(M)!=0]));
}

#' Network SD Vulnerability
#' @template M
SDVulnerability <- function(M){
  return(sd(OutDegree(M)[OutDegree(M)!=0]));
}

#' Network normalised mean Generality
#' @template M
MeanGenerality_norm <- function(M){
  norm_g <- NormalisedGenerality(M)
  return(mean(norm_g[norm_g!=0]))
}

#' Network normalised mean Generality -- following Williams and Martinez 2000
#' @template M
MeanGenerality_norm2 <- function(M){
  norm_g <- NormalisedGenerality(M)   ##  follows Williams and Martinez 2000
  return(mean(norm_g))
}

#' Network normalised mean Vulnerability
#' @template M
MeanVulnerability_norm <- function(M){
  norm_v <- NormalisedVulnerability(M)
  return(mean(norm_v[norm_v!=0]))
}

#' Network normalised mean Vulnerability -- following Williams and Martinez 2000
#' @template M
MeanVulnerability_norm2 <- function(M){
  norm_v <- NormalisedVulnerability(M)   ##  follows Williams and Martinez 2000
  return(mean(norm_v))
}

#' Network normalised SD Generality
#' @template M
SDGenerality_norm <- function(M){
  norm_g <- NormalisedGenerality(M)
  return(sd(norm_g[norm_g!=0])) ## this doesn't follow Williams and Martinez 2000 as it should
}

#' Network normalised SD Generality -- following Williams and Martinez 2000
#' @template M
SDGenerality_norm2 <- function(M){
  norm_g <- NormalisedGenerality(M)
  return(sd(norm_g)) ##  follows Williams and Martinez 2000
}

#' Network normalised SD Vulnerability
#' @template M
SDVulnerability_norm <- function(M){
  norm_v <- NormalisedVulnerability(M)
  return(sd(norm_v[norm_v!=0]))   ## this doesn't follow Williams and Martinez 2000 as it should
}

#' Network normalised SD Vulnerability -- following Williams and Martinez 2000
#' @template M
SDVulnerability_norm2 <- function(M){
  norm_v <- NormalisedVulnerability(M)
  return(sd(norm_v))   ##  follows Williams and Martinez 2000
}

#' Proportion of basal nodes
#' @template M
FractionOfBasal <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;

  b_sps <- sum(which(InDegree(M_temp) == 0) %in% which(OutDegree(M_temp) >= 1));

  return(b_sps / dim(M)[1]);
}

#' Number of basal nodes
#' @template M
NumberOfBasal <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;

  b_sps <- sum(which(InDegree(M_temp) == 0) %in% which(OutDegree(M_temp) >= 1));

  return(b_sps);
}

#' Proportion of top consumers
#' @template M
FractionOfTop <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;

  t_sps <- sum(which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) == 0));

  return(t_sps / dim(M)[1]);
}

#' Number of top consumers
#' @template M
NumberOfTop <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;

  t_sps <- sum(which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) == 0));

  return(t_sps);
}

# indegree_top <- function(M){
#   M_temp <- M;
#   diag(M_temp) <- 0;
#
#   top_indegree <- sum(M_temp[which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) == 0)]);
#   t_sps <- sum(which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) == 0));
#   t_indegree <- top_indegree/t_sps
#
#   return(t_indegree)
# }

#' Proportion of intermediate consumers
#' @template M
FractionOfIntermediate <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;

  i_sps <- sum(which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) >= 1));

  return(i_sps / dim(M)[1]);
}

#' Number of intermediate consumers
#' @template M
NumberOfIntermediate <- function(M){
  M_temp <- M;
  diag(M_temp) <- 0;

  i_sps <- sum(which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) >= 1));

  return(i_sps);
}

# indegree_intermediate <- function(M){
#   M_temp <- M;
#   diag(M_temp) <- 0;
#
#   intermediate_indegree <- sum(M_temp[which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) >= 1)]);
#   i_sps <- sum(which(InDegree(M_temp) >= 1) %in% which(OutDegree(M_temp) >= 1));
#   i_indegree <- intermediate_indegree/i_sps
#
#   return(i_indegree)
# }

#' Proportion of cannibalistic links
#' @template M
FractionOfCannibalism <- function(M){
  return(sum(diag(M)) / dim(M)[1]);
}

#' Network similarity
#' this is more like mean similarity than "maximum"
#' @template M
MaximumSimilarity <- function(M){
  S <- dim(M)[1];

  similarity <- 0;

  for(i in 1:S){
    for(j in 1:S){
      if(i == j) next;

      similarity <- similarity + ((sum(M[,i] & M[,j]) + sum(M[i,] & M[j,])) / (sum(M[,i] | M[,j]) + sum(M[i,] | M[j,])));

    }

  }
  return(similarity/S);
}

# From Nuria, but cheddar omnivory function doesn't deal well with loops
# Omnivory <- function(M){
#   #Use PredationMatrixToLinks() to create a Cheddar community from a predation
#   community <- Community(data.frame(node = colnames(M)), trophic.links = PredationMatrixToLinks(M),
#                          properties = list(title = "Test2"))
#   community <- RemoveCannibalisticLinks(community, title='community');
#
#   #community is a cheddar community
#
#   Fractionomnivory<-FractionOmnivorous(community)
#
#   return(Fractionomnivory)
# }

#' Network omnivory
#' @template M
#' @template dietcat
#' @param level a function to calculate the trophic level. See \code{cheddar}
#'   for options
Omnivory2 <- function(M, dietcat = NULL, level = PreyAveragedTrophicLevel){
  omnivs <- IsOmnivore2(M, level = level)
  if(!is.null(dietcat)){
    omnivs <- omnivs[!names(omnivs) %in% dietcat]
  }
  return(sum(omnivs) / length(omnivs))

}

#' Network mean food chain length
#' @template M
#'
#' @importFrom cheddar TrophicChainsStats Community RemoveCannibalisticLinks PredationMatrixToLinks
MeanFoodChainLength <- function(M){

  # Use PredationMatrixToLinks() to create a Cheddar community from a predation
  # matrix

  community <- Community(data.frame(node = colnames(M)), trophic.links = PredationMatrixToLinks(M),
                         properties = list(title = "Test2"))
  community <- RemoveCannibalisticLinks(community, title='community');

  # community is a Cheddar community
  #community
  #NPS(community)
  #TLPS(community)
  #TrophicLevels(community)

  #You can add node properties such as category:
  #category <- c('producer', 'invertebrate', 'vert.endo')
  #community <- Community(nodes=data.frame(node=node, category=category),
  #                       trophic.links=PredationMatrixToLinks(pm),
  #                       properties=list(title='Test community'))
  #NPS(community)

  #chs <- TrophicChains(community);
  #ch_lens <- ChainLength(chs);

  chain.stats <- TrophicChainsStats(community)
  ch_lens <- (chain.stats$chain.lengths + 1)

  return(sum(ch_lens)/length(ch_lens));
}

#' Network average predator overlap
#' @template M
CalculatePredatorOverlap <- function(M){
  cols <- dim(M)[2];
  M1 <- matrix(0, cols, cols);
  for(r in 1:dim(M)[1]){
    links <- which(M[r,] == 1);
    for(l in links){
      M1[l, links] <- 1;
    }
  }

  M1[lower.tri(M1)] <- 0;
  diag(M1) <- 0;

  return((sum(M1)) / ((cols) * ((cols-1)/2)));
}


################## END NETWORK METRICS ##################

#' Density of degree distribution calculation
#'
#' @template M
#' @param cumulative cumulative degree distribution?
#' @param ... passed to \code{igraph::degree}
#'
#' @importFrom graphics hist
#' @importFrom igraph is.igraph degree
#' @export
densityDegreeDistribution <- function(M, cumulative = TRUE, ...){
  graph <- graph.adjacency(M)
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  cs <- degree(graph, ...)
  hi <- hist(cs, -1:max(cs), plot = FALSE)$density
  if (!cumulative) {
    res <- hi
  }
  else {
    res <- rev(cumsum(rev(hi)))
  }
  res
}

