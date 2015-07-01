library(igraph)
library(reshape2)
library(data.table)
library(ggplot2)
library(Rcpp)
library(rbenchmark)
library(scales)
library(pii)
library(dplyr)
library(magrittr)
library(gridExtra)
library(parallel)

options(stringsAsFactors=F)

randomGraph <- function(i=0){
  graph <- watts.strogatz.game(1, 30, 3, 0.05)
  n <- sample(1:20, 1) / 100
  E(graph)$valence <- sample(c(-1, 1), ecount(graph), replace = T, prob = c(n, 1-n))
  E(graph)$color <- ifelse(E(graph)$valence == -1, "red", "black")
  graphid <- paste0('watts-strogatz ', '#',i, sep='')
  graph <- set.graph.attribute(graph, 'graphid', graphid)
  graph
}
g <- randomGraph()

getAllRingGraphs <- function(num = 5, chain = F){
  g <- graph.ring(num)
  if(chain) g <- delete.edges(g, 1)
  vals <- list()
  length(vals) <- length(E(g))
  for(i in 1:length(vals)){
    vals[[i]] <- c(1, -1)
  }
  all.combos <- expand.grid(vals)
  all.graphs <- list()
  for(r in 1:length(all.combos[,1])){
    E(g)$valence <- as.integer(all.combos[r,])
    E(g)$color <- ifelse(E(g)$valence == -1, "red", "black")
    V(g)$name <- letters[1:vcount(g)]
    graphid <- paste0(ifelse(chain, 'chain',  'ring'),
                      num, 'valence', paste0(as.integer(all.combos[r,]), collapse=''))
    g <- set.graph.attribute(g, 'graphid', graphid)
    all.graphs[[graphid]] <- g
  }
  return(all.graphs)
}
all.ring.graphs <- getAllRingGraphs()


graphData <- function(g) {
  data.table(graphid = get.graph.attribute(g, 'graphid'),
             meanDegree = mean(degree(g)),
             density = graph.density(g),
             diameter = diameter(g),
             nodeCount = vcount(g),
             edgeCount = ecount(g),
             degCentralization = centralization.degree(g)$centralization,
             avgPathLength = average.path.length(g),
             numPosEdge = length(E(g)[valence == 1]),
             numNegEdge = length(E(g)[valence == -1]),
             propNegEdge = length(E(g)[valence == -1]) / ecount(g),
             modularity = modularity(walktrap.community(g)),
             rankCor = suppressWarnings(cor(pii(g, pii.beta = -0.9), pii(g, pii.beta=-0.5), method = "spearman")),
             meanTrans = mean(transitivity(g, type='local'), na.rm=T))
}

nodeData <- function(g){
  beta.sequence <- seq(-0.9, -0.5, by=0.05)
  all.pii <- data.table(node = numeric(), pii.value = numeric(), degree = numeric(), beta = numeric(), graphid = character(), nodeid = character())
  g.ed <- edge.distance(g)
  g.gid <- get.graph.attribute(g, 'graphid')
  for(b in beta.sequence) {
    g.pii <- pii(g, e.dist = g.ed, pii.beta = b)
    td <- data.table(node=1:vcount(g), pii.value=as.numeric(g.pii), degree = degree(g), beta=b, graphid = g.gid, nodeid = letters[1:vcount(g)])
    all.pii <- rbind(all.pii, td)
  }
  return(all.pii)
}

lowCor <- function(g){
  p1 <- pii(g, pii.beta = -1)
  breakval = 0
  for(b in seq(-0.99, -0.01, by=0.01)){
    p <- pii(g, pii.beta = b)
    rc <- cor(p1, p, method = "spearman")
    if(rc <= 0.707){
      breakval = b
      break
    }
  }
  return(breakval)
}

avgMinDistNegEdge <- function(g){
  negEdgeDist <- edge.distance(g)[which(E(g)$valence == -1), ]
  x <- 0
  for(i in 1:length(negEdgeDist[1,])){
    x <- x + min(negEdgeDist[,i])
  }
  return(x / length(negEdgeDist))
}

avgDistNegEdge <- function(g){
  negEdgeDist <- edge.distance(g)[which(E(g)$valence == -1), ]
  x <- 0
  for(i in 1:length(negEdgeDist[1,])){
    x <- x + sum(negEdgeDist[,i])
  }
  return(x / length(negEdgeDist))
}

