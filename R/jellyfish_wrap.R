
jelly_count <- function(x, k = 20, hash = 10^8, threads = 1, jellyfish.exec = "jellyfish"){
  
  if(is(x, "DNAbin")){
    fa_file <- tempfile()
    write.dna(x, file = fa_file, format = "fasta")
  } else if(file.exists(x)){
    fa_file <- normalizePath(x)
    fa_file <- gsub(" ", "\\ ", fa_file, fixed = TRUE)
  } else {
    stop("Invalid input: x must be a DNAbin object or a file path.")
  }
  
  jf_file <- tempfile()
  ct_file <- tempfile()
  hash <- format(hash, scientific = FALSE)
  
  comm <- paste(jellyfish.exec, "count",
                "-m", k,
                "-t", threads,
                "-s", hash,
                "-C", fa_file,
                "-o", jf_file, "&&",
                jellyfish.exec, "dump",
                jf_file, ">", ct_file)
  system(comm)
  
  res <- read.FASTA(ct_file)
  res <- as.character(res)
  res <- sapply(res, function(x){paste(toupper(x), collapse = "")})
  res <- data.frame(kmer = res, n = as.numeric(names(res)), stringsAsFactors = FALSE)
  
  file.remove(jf_file, ct_file)
  if(is(x, "DNAbin")){
    file.remove(fa_file)
  }
  
  return(res)
}




