#' diatobc
#' @importFrom ape read.dna write.dna as.DNAbin labels.DNAbin as.matrix.DNAbin rbind.DNAbin read.FASTA read.tree write.tree multi2di drop.tip is.binary.tree muscle 
#' @importFrom phylobase tipLabels
#' @importFrom adephylo distTips
#' @importFrom jsonlite fromJSON
#' @importFrom ips phylo2mafft write.fas read.fas
#' @importFrom stringr str_detect
#' @importFrom seqinr consensus
#' @importFrom methods is
#' 
#' 
#' 
#' @importFrom grDevices topo.colors
#' @importFrom graphics plot
#' @importFrom stats quantile
#' @importFrom utils setTxtProgressBar txtProgressBar
#' 
#' @name diatobc
#' @docType package
NULL



#' Diatoms Phylogeny
#' 
#' A phylogenetic tree of 548 diatom species.
#' 
#' @format a phylogenetic tree in \code{phylo} format.
#' @docType data
#' @keywords datasets
#' @name diatoms.tree
#' @usage data(diatoms.tree)
#' @details The phylogeny has been reconstructed with RaxML and two genetic markers (18S and rbcl).
#' The tree is not rooted but can be rooted with 'Bolidomonas pacifica', a non-diatom species included in the phylogeny.
NULL



#' Diatoms genetic sequences
#' 
#' An alignment of genetic sequences for two markers (18S and rbcl) for
#' 548 diatom species and a non-diatom species ('Bolidomonas pacifica').
#' 
#' @format an alignment in \code{DNAbin} format.
#' @docType data
#' @keywords datasets
#' @name diatoms.seq
#' @usage data(diatoms.seq)
NULL


#' Diatoms IPS-P68 index
#' 
#' This dataset contains all the necessary data to compute the IPS-P index based on 68 clusters as described in Keck et al. (2016).
#' 
#' @format a list of two elements:
#' \describe{
#'   \item{clusters}{a \code{graphclust} object (package \pkg{phylosignal}) giving the 68 clusters of the index.}
#'   \item{scores}{IPSS (sensitivity) and IPSV (variability) for each cluster}
#' }
#' @references Keck F., Rimet F., Franc A. & Bouchez A. (2016) Linking Diatoms Ecological Preferences To Phylogeny: New Perspectives for Aquatic Ecosystems Bioassessment.
#' @docType data
#' @keywords datasets
#' @name IPSP68
#' @usage data(IPSP68)
NULL