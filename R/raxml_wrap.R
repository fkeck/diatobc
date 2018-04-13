
#' RAxML phylogeny with backbone tree. 
#' 
#' This function uses RAxML to reconstruct a tree given an alignment
#' and a backbone tree (the -g option of RAxML).
#' 
#' @param x an alignment of sequences of class \code{DNAbin}.
#' @param backbone a phylogenetic tree (class \code{phylo}) to use as a constraint.
#' The tip labels must match with the appropriate DNA sequences labels. Use \code{NULL} for no constraint.
#' @param model a character string giving the substitution model to use (the -m option of RAxML). 
#' @param raxml.exec a character string giving the path of the program.
#' @param nruns the number of independent ML searches (runs). This is ignored if bootstrapping is enabled.
#' Then the number of ML searches will be the number of bootstraps set with \code{bs} divided by 5.
#' @param bs a numeric giving the number of rapid bootstraps to perform.
#' @param part a vector of integers giving the boundaries of the different partitions, if relevant.
#' E.g. \code{c(1, 5355, 5356, 7633)} for a two genes partition.
#' @param threads the number of threads to use, if available.
#' @param MLSCC a logical indicating whether to use the ML search convergence criterion (the -D option of RAxML).
#' This is a shortcut in ML search which can improve significantly the speed for large phylogenies.
#' @param quiet a logical indicating whether RAxML warnings are printed on console (available only for RAxML v8+).
#' 
#' @details Since RAxML doesn't support special characters in taxa/sequence names,
#' they are deleted automatically with a warning.
#' Be warned that deletion of illegal characters may lead to ambiguities in taxa/sequence names.
#' 
#' @return A \code{list} of two elements:
#' \describe{
#'   \item{\code{tree}}{The best-scoring ML tree as a \code{phylo} object}
#'   \item{\code{info}}{The log file printed by RAxML.}
#' }
#' 
#' @seealso \code{\link[ips]{raxml}} in the package \pkg{ips} for a general RAxML interface,
#' \code{\link{raxmlEPA}} for RAxML Evolutionary Placement Algorithm.
#' 
#' @references Stamatakis A. (2014) RAxML Version 8: A tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies. Bioinformatics.
#' @rdname raxml
#' @export
raxmlBackbone <- function(x, backbone, model, raxml.exec, nruns = 1, bs = NULL, part = NULL, threads = NULL, MLSCC = FALSE, quiet = FALSE){
  
  if(is.list(x)) x <- as.matrix(x)
  
  tmp <- paste0(tempdir(), "/raxmlBackbone")
  if(!file.exists(tmp)){
    dir.create(tmp)
  }
  file.remove(list.files(tmp, full.names = TRUE))
  
  input.seq <- paste0(tmp, "/input_seq.phy")
  input.tre <- paste0(tmp, "/input_tre.tre")
  input.par <- paste0(tmp, "/input_par.txt")
  if(is.null(bs)){
    output.tre <- paste0(tmp, "/RAxML_bestTree.diatobc_out")
  } else {
    output.tre <- paste0(tmp, "/RAxML_bipartitions.diatobc_out")
  }
  output.inf <- paste0(tmp, "/RAxML_info.diatobc_out")
  
  valid.labs <- rownames(x)
  names(valid.labs) <- gsub("[^0-9a-zA-Z_]+", "", rownames(x))
  
  rownames(x) <- gsub("[^0-9a-zA-Z_]+", "", rownames(x))
  write.dna(x, file = input.seq, format = "interleaved")
  
  if(is.null(bs)){
    com.mode <- "-f d"
  } else {
    com.mode <- "-f a -x 1000"
    nruns <- bs
  }
  
  if(!is.null(backbone)){
    backbone$tip.label <- gsub("[^0-9a-zA-Z_]+", "", backbone$tip.label)
    write.tree(backbone, file = input.tre)
    com.backbone <- paste("-g", input.tre)
  } else {
    com.backbone <- ""
  }

  
  if(!is.null(part)){
    part <- raxmlPart(part)
    write(part, file = input.par)
    com.part <- paste("-q", input.par)
  } else {
    com.part <- ""
  }
  
  if(!is.null(threads)){
    com.threads <- paste("-T", threads)
  } else {
    com.threads <- ""
  }
  
  if(MLSCC){
    com.mlscc <- "-D"
  } else {
    com.mlscc <- ""
  }
  
  if(quiet){
    com.quiet <- "--silent"
  } else {
    com.quiet <- ""
  }
  
  
  com <- paste(raxml.exec,
               com.threads,
               com.mode,
               "-m", model,
               "-N", nruns,
               "-O",
               "-p 1000",
               "-s", input.seq,
               "-n diatobc_out",
               com.part,
               com.backbone,
               "-w", tmp,
               com.mlscc,
               com.quiet)
  
  com <- gsub(" +", " ", com)
  
  system(com)
  
  res <- list()
  res$tree <- read.tree(output.tre)
  res$tree$tip.label <- valid.labs[res$tree$tip.label]
  res$info <- scan(file = output.inf, what = "character", sep = "\n", quiet = TRUE)
  
  return(res)
}



#' @rdname raxml
#' @export
raxmlEPA <- function(x, ref.tree, model, raxml.exec, nruns = 1, part = NULL, threads = NULL, quiet = FALSE,
                     extra.com = ""){
  
  tmp <- paste0(tempdir(), "/raxmlEPA")
  if(!file.exists(tmp)){
    dir.create(tmp)
  }
  file.remove(list.files(tmp, full.names = TRUE))
  
  input.seq <- paste0(tmp, "/input_seq.phy")
  input.tre <- paste0(tmp, "/input_tre.tre")
  input.par <- paste0(tmp, "/input_par.txt")
  output.tre <- paste0(tmp, "/RAxML_labelledTree.diatobc_out")
  output.inf <- paste0(tmp, "/RAxML_info.diatobc_out")
  
  valid.labs <- rownames(x)
  names(valid.labs) <- gsub("[^0-9a-zA-Z_]+", "", rownames(x))
  rownames(x) <- gsub("[^0-9a-zA-Z_]+", "", rownames(x))
  write.dna(x, file = input.seq, format = "interleaved")
  
  com.mode <- "-f v"

  ref.tree$tip.label <- gsub("[^0-9a-zA-Z_]+", "", ref.tree$tip.label)
  write.tree(ref.tree, file = input.tre)
  com.ref.tree <- paste("-t", input.tre)
  
  if(!is.null(part)){
    part <- raxmlPart(part)
    write(part, file = input.par)
    com.part <- paste("-q", input.par)
  } else {
    com.part <- ""
  }
  
  if(!is.null(threads)){
    com.threads <- paste("-T", threads)
  } else {
    com.threads <- ""
  }
  
  if(quiet){
    com.quiet <- "--silent"
  } else {
    com.quiet <- ""
  }
  
  com <- paste(raxml.exec,
               com.threads,
               com.mode,
               "-m", model,
               "-N", nruns,
               "-O",
               "-p 1000",
               "-s", input.seq,
               "-n diatobc_out",
               com.part,
               com.ref.tree,
               "-w", tmp,
               com.quiet,
               extra.com)
  
  com <- gsub(" +", " ", com)
  
  system(com)
  
  res <- list()
  res$tree <- read.tree(output.tre)
  res$tree$tip.label <- gsub("QUERY___", "", res$tree$tip.label)
  res$tree$tip.label <- valid.labs[res$tree$tip.label]
  res$info <- scan(file = output.inf, what = "character", sep = "\n", quiet = TRUE)
  
  return(res)
}
