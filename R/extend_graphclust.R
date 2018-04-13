#' Extend clusters to new species using a phylogenetic tree
#' 
#' This function affiliate new species to clusters of a \code{graphclust} object.
#' This is done using a phylogeny including both the species included in the \code{graphclust}
#' object and the set of new species to affiliate.
#' 
#' @param clusters an object of class \code{graphclust}
#' @param target a \code{phylo} object
#' @param pb a logical. Should a progress bar be printed? (default TRUE).
#' 
#' @details The affiliation algorithm is basic and works with phylogenetic distance as number of nodes between two tips.
#' It may produce \code{NA} values in case of ambiguity between two clusters.
#' Affliliation of a large number of new tips can be slow.
#' 
#' @return A vector giving the cluster of each tip.
#' 
#' @export
extendGraphClust <- function(clusters, target, pb = TRUE){
  
  tip.ref <- intersect(tipLabels(clusters$meta$p4d), target$tip.label)
  tip.assign <- setdiff(target$tip.label, tipLabels(clusters$meta$p4d))
  
  res <- vector(mode = "integer", length = length(tip.assign))
  
  if(pb){
    pb <- txtProgressBar(min = 0, max = length(tip.assign), style = 3)
  }
  
  for(i in 1:length(tip.assign)){
    tre.i <- drop.tip(target, tip.assign[-i])
    
    tre.i.dist <- as.matrix(distTips(tre.i, method = "nNodes"))
    tip.neig <- which(tre.i.dist[tip.assign[i],] <= min(tre.i.dist[tip.assign[i],][tre.i.dist[tip.assign[i],] > 0]) & tre.i.dist[tip.assign[i],] > 0)
    
    tip.neig.clu <- clusters$clusters[tip.neig]
    
    if(abs(max(tip.neig.clu) - min(tip.neig.clu)) == 0){
      res[i] <- tip.neig.clu[1]
    } else {
      res[i] <- NA
    }
    setTxtProgressBar(pb, i)
  }
  
  names(res) <- tip.assign
  return(res)
}