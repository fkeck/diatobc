#' Tree relative datation with PATHd8
#' 
#' This function uses PATHd8 to date phylogenetic trees with the Mean Path Length algorithm.
#' 
#' @param phy a phylogenetic tree (class \code{phylo}) to date.
#' @param pathd8.exec a character string giving the path of the program.
#' @details This is a limited wrapper of the C program PATHd8 to perform tree datation in relative time
#' with the Mean Path Length algorithm. Wrapper for calibrated datation (D8) is not implemented.
#' @return An ultrametric tree of class \code{phylo}.
#' @references
#' \itemize{
#'  \item Britton, T., Anderson, C. L., Jacquet, D., Lundqvist, S., & Bremer, K. (2007). Estimating divergence times in large phylogenetic trees. Systematic biology, 56(5), 741-752.
#'  \item Britton, T., Oxelman, B., Vinnersten, A., & Bremer, K. (2002). Phylogenetic dating with confidence intervals using mean path lengths. Molecular phylogenetics and evolution, 24(1), 58-65.
#' }
#' @export
#' 
pathd8 <- function(phy, pathd8.exec){
  tmp <- paste0(tempdir(), "/pathd8")
  if(!file.exists(tmp)){
    dir.create(tmp)
  }
  file.remove(list.files(tmp, full.names = TRUE))
  input.pd8 <- paste0(tmp, "/input_pd8.tre")
  output.pd8 <- paste0(tmp, "/output_pd8.txt")
  
  write.tree(phy = phy, file = input.pd8)
  com.pd8 <- paste(pathd8.exec, input.pd8, output.pd8)
  
  system(com.pd8)
  cat("Loading data...")
  res <- readLines(output.pd8)
  res <- res[grep("^d8 tree    : ", res)]
  res <- read.tree(text = res)
  cat("Done.")
  return(res)
}


