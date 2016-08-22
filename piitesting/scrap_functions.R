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
igraph.options(vertex.color = '#FFFFFFAA', vertex.alpha = 0.5, vertex.size = 15, vertex.frame.color = '#eeeeee',
               edge.width = 4, vertex.frame.size = 0,
               vertex.label.family = 'Futura Medium', vertex.label.cex = 1.5, vertex.label.color = '#333333')

thm <- theme(
  panel.background = element_rect(fill = 'white'),
  plot.background = element_rect(fill = 'white'),
  axis.ticks = element_blank(),
  axis.text = element_text(family = 'Futura Medium'),
  axis.title = element_text(family = 'Futura Medium', size = 15),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(color = '#AAAAAA')
)


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

guid <- function() {
  baseuuid <- paste(sample(c(letters[1:6],0:9),30,replace=TRUE),collapse="")

  paste(
    substr(baseuuid,1,8),
    "-",
    substr(baseuuid,9,12),
    "-",
    "4",
    substr(baseuuid,13,15),
    "-",
    sample(c("8","9","a","b"),1),
    substr(baseuuid,16,18),
    "-",
    substr(baseuuid,19,30),
    sep="",
    collapse=""
  )
}

randomGraph <- function(maxnegtie = 20, pendchance = 0.2, badeggchance = 0.05, benegpercent = 0.8){
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
      for(e in E(graph)[from(V(graph)[i])]) {
        E(graph)[e]$valence <- sample(c(-1, 1), 1, prob = c(benegpercent, 1 - benegpercent))
      }
    }
  }

  E(graph)$color <- ifelse(E(graph)$valence == -1, "red", "black")
  V(graph)$name <- letterNames(vcount(graph))
  graphid <- paste0('watts-strogatz ', '#', guid(), sep='')
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

# add number of crosses
# size-adjusted number of crosses
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
             optimBeta = optimBeta(g))
             #numCrosses = nrow(crosses(g)),
             #sizeAdjNumCross = nrow(crosses(g)) / factorial(vcount(g)))
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
  breakval = NA
  for(b in seq(-0.99, -0.01, by=0.01)){
    p <- pii(g, pii.beta = b)
    rc <- suppressWarnings(cor(p1, p, method = "spearman"))
    if(is.na(rc)) {
      return(NA)
    }
    if(rc <= 0.707) {
      return(b)
    }
  }
  return(breakval)
}

avgMinDistNegEdge <- function(g){
  if(clusters(g)$no > 1){return(NA)}          # if it's multiple components
  if(!any(E(g)$valence == -1)) { return(NA) } # if there's no neg edge
  negEdgeDist <- edge.distance(g)[which(E(g)$valence == -1), ]
  negEdgeDist <- as.matrix(negEdgeDist)
  x <- mean(apply(negEdgeDist, 2, min))
  return(x)
}

avgDistNegEdge <- function(g){
  if(clusters(g)$no > 1){return(NA)}
  if(!any(E(g)$valence == -1)) { return(NA) } # if there's no neg edge
  negEdgeDist <- edge.distance(g)[which(E(g)$valence == -1), ]
  negEdgeDist <- as.matrix(negEdgeDist)
  x <- sum(negEdgeDist)
  return(x / (vcount(g) * sum(E(g)$valence == -1)))
}

optimBeta <-  function(g, comp.left = -0.9, comp.right = -0.5, starting.beta = -0.8, unbounded = F, lower.bound = comp.left, upper.bound = comp.right) {
  rc.diff <- function(b, g, comp.left = comp.left, comp.right = comp.right) {
    p <- pii(g, pii.beta = b)
    pii.left <- pii(g, pii.beta = comp.left)
    pii.right <- pii(g, pii.beta = comp.right)
    rc.left = cor(p, pii.left, method = "spearman")
    rc.right = cor(p, pii.right, method = "spearman")
    return(abs(rc.left - rc.right))
  }

  if(unbounded) {
    o <- optim( par = c(starting.beta), fn = rc.diff, gr = NULL,
                g = g, comp.right = comp.right, comp.left = comp.left, method = 'Brent',
                lower = -1, upper = -0.0001)
  } else {
    o <- optim( par = c(starting.beta), fn = rc.diff, gr = NULL,
                g = g, comp.right = comp.right, comp.left = comp.left, method = 'Brent',
                lower = lower.bound, upper = upper.bound)
  }

  cross.point <- o$par
  p <- pii(g, pii.beta = o$par)
  pii.left <- pii(g, pii.beta = comp.left)
  attr(cross.point, 'spearman.correlation') <- cor(p, pii.left, method = "spearman")
  attr(cross.point, 'approx') <- o$value
  return(cross.point)
}

pii.diagnostic.plot <- function(g) {
  require(ggplot2)
  require(parallel)
  require(gridExtra)

  inc <- 0.01
  x <- do.call('rbind', mclapply(seq(-1, -0.01, by=inc), function(b) {
    p <- pii(g, pii.beta = b)
    piir <- rank(p)
    data.table(b = b, nd = V(g)$name, pii = p, piir = piir)
  }, mc.cores = 6))

  p1 <- ggplot(x, aes(x=b, y=pii, group=nd, color=nd)) +
    scale_x_continuous(expression(beta), limits = c(-1,-0.01)) +
    scale_y_continuous('PII Score') +
    geom_line(size=1) + thm + theme(legend.position = 'none')

  # color the lines by pii rank change
  piic <- x %>% filter(b %in% c(-0.9, -0.5)) %>%
    dcast(nd ~ b, value.var = 'piir') %>%
    mutate(piic = `-0.9` - `-0.5`)
  x <- x %>% left_join(piic, by = 'nd')
  p2 <- ggplot(x, aes(x=b, y=piir, group=nd, color=piic)) +
    geom_line(alpha=0.9, size = 1) +
    scale_color_distiller(type = 'div', palette = 5) +
    scale_x_continuous(expression(beta), lim = c(-1, -0.01)) +
    scale_y_continuous('PII Ranked Score') +
    thm + theme(legend.position = 'none', panel.grid.major.y = element_blank())

  comp.left <- -0.9
  comp.right <- -0.5
  rc <- do.call('rbind', lapply(unique(x$b)[-1], function(bb) {
    pii <- x[b == bb, pii]
    pii.left <- x[b == comp.left, pii]
    pii.right <- x[b == comp.right, pii]
    data.table(b = bb, rc.left = cor(pii, pii.left, method = "spearman"), rc.right = cor(pii, pii.right, method = "spearman"))
  }))

  # rc.diff <- function(b, g, comp.left=-0.9, comp.right=-0.5) {
  #   p <- pii(g, pii.beta = b)
  #   pii.left <- pii(g, pii.beta = comp.left)
  #   pii.right <- pii(g, pii.beta = comp.right)
  #   rc.left = cor(p, pii.left, method = "spearman")
  #   rc.right = cor(p, pii.right, method = "spearman")
  #   return(abs(rc.left - rc.right))
  # }
  # rc.diff(-0.63, g)
  # o <- optim(par=c(-0.8), fn = rc.diff, gr = NULL,
  #            g = g, comp.right = comp.right, comp.left = comp.left,
  #            method='Brent', lower = -1, upper = -0.001)
  # cross.y <- cor(pii(g, pii.beta = o$par), pii(g, pii.beta = comp.left), method = "spearman")
  # cross.point <- data.frame(xx=o$par, yy=cross.y)

  ## Need to find the crossing
  rc.m <- melt(rc, id.vars='b')
  p3 <- ggplot() +
    geom_line(data = rc.m, aes(x=b, y=value, group=variable)) +
    # geom_text(data=cross.point, aes(x=xx, y=yy, label = round(xx,3)), size = 4, vjust=2) +
    scale_y_continuous('Rank Correlation', lim=c(0,1.0)) +
    scale_x_continuous(expression(beta), lim=c(-1, -0.01)) + thm

  grid.arrange(p2, p3)
}

