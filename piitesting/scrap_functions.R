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
library(parallel)

options(stringsAsFactors=F)
igraph.options(vertex.color = '#FFFFFFAA', vertex.alpha = 0.5, vertex.size = 15, vertex.frame.color = 'grey60',
               vertex.label.family = 'Ubuntu', vertex.label.cex = 0.7, vertex.label.color = 'blue')



letterNames <- function(n) {
  nms <- character(n)
  k <- 1
  for(i in 1:26) {
    for(j in 1:26) {
      nms[k] <- paste0(letters[i], letters[j])
      k <- k + 1
      if(k > n) {
        return(nms)
      }
    }
  }
}


# change the percentages to actual percentages using runif(1) < p

randomGraph <- function(i=0, maxnegtie = 20, pendchance = 0.2, badeggchance = 0.05, benegpercent = 80){
  N <- sample(10:100, 1)
  graph <- watts.strogatz.game(1, N, 3, 0.05)

  #random pendants
  for(i in 1:length(V(graph))){
    if(runif(1) < pendchance){
      graph <- add.vertices(graph, 1)
      graph <- add.edges(graph, c(V(graph)[i], V(graph)[length(V(graph))]))
    }
  }
  n <- sample(1:maxnegtie, 1) / 100

  E(graph)$valence <- sample(c(-1, 1), ecount(graph), replace = T, prob = c(n, 1-n))

  #bad eggs
  for(i in 1:length(V(graph))){
    if(runif(1) < badeggchance){
      E(graph)[from(V(graph)[i])]$valence <- -1
    }
  }

  E(graph)$color <- ifelse(E(graph)$valence == -1, "red", "black")
  V(graph)$name <- letterNames(vcount(graph))
  graphid <- paste0('watts-strogatz ', '#',i, sep='')
  graph <- set.graph.attribute(graph, 'graphid', graphid)
  graph
}

modularGraph <- function(nummods = 2, orignummods = nummods, numconedges = 3, graph = NULL){
  g1 <- randomGraph()
  V(g1)$name <- paste0(V(g1)$name, '1')
  V(g1)$grp <- 1
  if(!(is.null(graph))){
    g1 <- graph
  }
  #set.vertex.attribute(g1, name = 'graphname', value = '1')
  g2 <- randomGraph()
  numdiff <- orignummods - nummods + 2
  V(g2)$name <- paste0(V(g2)$name, numdiff)
  V(g2)$grp <- numdiff
  #set.vertex.attribute(g1, name = 'graphname', value = '2')
  igraph.options(edge.attr.comb = list(valence = max))

  e1 <- get.data.frame(g1)
  v1 <- get.data.frame(g1, what = 'vertices')

  e2 <- get.data.frame(g2)
  v2 <- get.data.frame(g2, what = 'vertices')

  e.all <- rbind(e1, e2)
  v.all <- rbind(v1, v2)

  g3 <- graph.data.frame(e.all, directed = F, vertices = v.all)

  if(nummods > 2){
    nummods_2 <- nummods-1
    g3 <- modularGraph(nummods = nummods_2, orignummods = orignummods, graph = g3)
  }

  return(g3)
  #g3 <- graph.union(g1, g2)

  #rand1 <- sample(V(g3)[V(g3)$grp_1 == 1], 1)
  #rand2 <- sample(length(V(g1):length(V(g2)), 1)

  #g3 <- add.edges(g3, c(V(g3)[rand1], V(g3)[rand2]))
  #E(g3)[length(E(g3)[from(rand1)])]$valence <- 1

  #E(g3)$color <- ifelse(!is.na(E(g3)$color_1), E(g3)$color_1, E(g3)$color_2)
  #for(i in 1:length(E(g3))){
  #  if(!is.na(E(g3)[i]$valence)){
  #    E(g3)[i]$color <- ifelse(E(g3)[i]$valence == -1, "red", "black")
  #  }
  #}


}

tieAdder <- function(g, tiechance = 0.001){
  e <- get.data.frame(g)
  v <- get.data.frame(g, what = 'vertices')

  for(i in 1:length(v[,1])){
    for(j in (i+1):length(v[,1])){
      if(v$grp[j] != v$grp[i] & runif(1) < tiechance){
        val <- sample(c(1,-1), 1, prob = c(0.2, 0.8))
        col <- ifelse(val == -1, "red", "black")

        e <- rbind(e,
                 data.frame(from = v$name[i], to = v$name[j],
                            valence=val, color = col))
      }
    }
  }
  g <- graph.data.frame(e, directed = F, vertices = v)
  return(g)
}

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
             meanTrans = mean(transitivity(g, type='local'), na.rm=T),
             lowCorBeta = lowCor(g),
             avgMinDistToNegEdge = avgMinDistNegEdge(g),
             avdDistOfNegEdge = avgDistNegEdge(g),
             optimBeta(g))

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
  if(clusters(g)$no > 1){return(NA)}
  negEdgeDist <- edge.distance(g)[which(E(g)$valence == -1), ]
  x <- 0
  for(i in 1:ncol(negEdgeDist)) {
    x <- x + min(negEdgeDist[,i])
  }
  return(x / vcount(g))
}

avgDistNegEdge <- function(g){
  if(clusters(g)$no > 1){return(NA)}
  negEdgeDist <- edge.distance(g)[which(E(g)$valence == -1), ]
  x <- 0
  for(i in 1:ncol(negEdgeDist)) {
    x <- x + sum(negEdgeDist[,i])
  }
  return(x / vcount(g))
}

optimBeta <- function(g, init.beta = -0.8, comp.left=-0.9, comp.right=-0.5) {

  rc.diff <- function(b, g) {
    p <- pii(g, pii.beta = b)
    pii.left <- pii(g, pii.beta = comp.left)
    pii.right <- pii(g, pii.beta = comp.right)
    rc.left = cor(p, pii.left, method = "spearman")
    rc.right = cor(p, pii.right, method = "spearman")
    return(abs(rc.left - rc.right))
  }

  o <- optim(par= init.beta, fn = rc.diff, gr = NULL, g = g,
             method='Brent', lower = -1, upper = -0.001)
  return(o$par)
}

