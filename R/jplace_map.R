






jplaceMap <- function(input){
  
  options(scipen=999)
  json <- jsonlite::fromJSON(input, simplifyDataFrame = FALSE)
  
  tree <- json$tree
  tree <- gsub("\\{", "\\[", tree)
  tree <- gsub("\\}", "\\]", tree)
  trr <- read.tree(text = tree)

  tree2 <- gsub(":(.*?)\\[", ":[", tree)
  tree2 <- gsub("\\]|\\[", "", tree2)
  trr2 <- read.tree(text = tree2)
  trr2$edge.length <- trr2$edge.length + 1
  
  pqueries.values <- lapply(json$placements, function(x) x[["p"]])
  pqueries.values <- do.call(rbind, pqueries.values)
  colnames(pqueries.values) <- json$fields
  pqueries.values[, "edge_num"] <- pqueries.values[, "edge_num"] + 1
  
  pqueries.values[, "edge_num"] <- match(pqueries.values[, "edge_num"], trr2$edge.length)
  
  dl.edge <- tapply(pqueries.values[ , "distal_length"], pqueries.values[ , "edge_num"], c)
  
  trr.c <- pathd8(trr, pathd8.exec = "/home/regalec/Logiciels/PATHd8/PATHd8")
  trr.c$edge.length <- trr.c$edge.length + 0.000003
  ratio.ultra <- (trr.c$edge.length / trr$edge.length)[as.numeric(names(dl.edge))]
  dl.edge <- mapply(function(x, y) x * y, x = dl.edge, y = ratio.ultra)
  
  resolution <- 100
  el <- trr.c$edge.length - 0.000002
  fragment.unit <- max(trr.c$edge.length)/resolution
  maps <- lapply(el, function(x) c(0.000001, rep(fragment.unit, x %/% fragment.unit), x %% fragment.unit, 0.000001))
  
  count.insert <- function(x){
    maps.data <- rep(0, length(maps[[x]]))
    counted <- table(findInterval(dl.edge[[as.character(x)]], cumsum(maps[[x]]), all.inside = TRUE))
    if(length(maps.data) == 1){
      counted <- sum(counted)
      names(counted) <- 1
    }
    maps.data[as.numeric(names(counted)) + 1] <- counted
    return(maps.data)
  }
  
  maps.data <- lapply(1:length(trr.c$edge.length), count.insert)
  maps.data <- lapply(maps.data, function(x) log(x+1.1))
  coef.ncol <- 1000 / max(unlist(maps.data))
  maps.data <- lapply(maps.data, function(x) round(x * coef.ncol, digits = 0))
  maps.data <- lapply(maps.data, function(x) {y <- x ; y[x == 0] <- 1; return(y)})
  
  maps <- mapply(function(x, y) {z <- x; names(z) <- as.character(y); return(z)}, x = maps, y = maps.data)
  
  
  cols <- topo.colors(1000)
  names(cols) <- 1:1000
  res <- list(tree = trr.c, cols = cols, lims = c(0, 1))
  class(res) <- "contMap"
  res$tree$maps <- maps
  
  res$tree$tip.label <- gsub("_", " ", res$tree$tip.label)
  #res$tree$tip.label <- binom2omnidia(res$tree$tip.label)
  return(res)
}


#input <- "Results_BC_insertion/EPA/RAxML_portableTree.diatobc_out.jplace"

#phytools::plot.contMap(res, type = "fan", fsize = 0.5, legend = FALSE)

