#' Insert barcodes in an alignment
#' 
#' This function inserts new barcodes into an alignment.
#' For this purpose, the barcode region of the reference alignment is re-aligned with the new set of barcodes.
#' 
#' @param bc barcodes sequences. A \code{DNAbin} object or a path of a sequence file to be read in fasta format.
#' @param ref the reference alignment. A \code{DNAbin} object or a path of a sequence file to be read in fasta format.
#' @param region a vector of two integers giving the boundaries of the barcodes region within the reference alignment.
#' @param align.method a character string giving the method to use for alignment.
#' Available methods are "\code{mafft}", "\code{muscle}" and "\code{none}" (to avoid alignment step).
#' @param align.exec a character string giving the path of the program if necessary.
#' @param threads an integer giving the number of threads to use. Available only with mafft.
#' @param trim logical. If \code{TRUE}, only the barcode region is returned.
#' 
#' @return A \code{DNAbin} alignment.
#' 
#' @export
insertBC <- function(bc, ref, region, align.method = "mafft", align.exec = NULL, threads = NULL, trim = FALSE){
  
  align.method <- match.arg(align.method, c("mafft", "muscle", "none"), several.ok = FALSE)
  if(is.null(align.exec)){
    align.exec <- align.method
  }
  if(!is(ref, "DNAbin")){
    ref <- read.dna(ref, format = "fasta")
  }
  if(!is.matrix(ref)){
    stop("ref sequences must be aligned")
  }
  
  if(!is(bc, "DNAbin")){
    bc <- read.dna(bc, format = "fasta")
  }
  bc <- as.list(bc)
  
  ref.a <- ref[, 1:(min(region)-1)]
  ref.b <- ref[, min(region):max(region)]
  ref.c <- ref[, (max(region)+1):ncol(ref)]
  
  ref.b.bc <- c(as.list(ref.b), bc)
  if(align.method == "mafft"){
    if(!is.null(threads)){
      threads <- paste("--thread", threads)
    } else {
      threads <- ""
    }
    ref.b.bc.align <- mafft2(x = ref.b.bc, path = align.exec, quiet = FALSE, options = threads)
  }
  
  if(align.method == "muscle"){
    ref.b.bc.align <- muscle(ref.b.bc, exec = align.exec, quiet = FALSE, original.ordering = TRUE)
  }
  
  if(align.method == "none"){
    ref.b.bc.align <- as.matrix(ref.b.bc)
  }
  
  labels.empty <- setdiff(labels(ref.b.bc), labels(ref.b.bc.align))
  ref.b.bc.empty <- matrix("-", nrow = length(labels.empty), ncol = ncol(ref.b.bc.align), dimnames = list(labels.empty))
  ref.b.bc.empty <- as.DNAbin(ref.b.bc.empty)
  ref.b.bc.align <- rbind(ref.b.bc.align, ref.b.bc.empty)
  ref.b.bc.align <- ref.b.bc.align[names(ref.b.bc), ]
  
  ref.a.bc.empty <- matrix("-", nrow = length(bc), ncol = ncol(ref.a), dimnames = list(labels(bc)))
  ref.a.bc.empty <- as.DNAbin(ref.a.bc.empty)
  ref.a.bc.align <- rbind(ref.a, ref.a.bc.empty)
  
  ref.c.bc.empty <- matrix("-", nrow = length(bc), ncol = ncol(ref.c), dimnames = list(labels(bc)))
  ref.c.bc.empty <- as.DNAbin(ref.c.bc.empty)
  ref.c.bc.align <- rbind(ref.c, ref.c.bc.empty)
  
  if(trim){
    res <- ref.b.bc.align
  } else {
    res <- cbind(ref.a.bc.align, ref.b.bc.align, ref.c.bc.align, check.names = FALSE)
  }
  return(res)
}
