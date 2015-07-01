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

options(stringsAsFactors=F)

###### trim ##############################
# trims the leading and trailing spaces from a string (you can provide a vector)
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

plotty <- function(g) {
  if(!('igraph' %in% class(g))) {
    stop('Requires an igraph object.')
  } else {
    igraph::plot.igraph(g, vertex.size=8, vertex.color='SkyBlue2', vertex.frame.color='SkyBlue2', edge.arrow.size=0.5,
         vertex.label.family='Open Sans', vertex.label.cex=0.5)
  }
}

# create some random graph
N <- 10
g <- erdos.renyi.game(N, runif(N))

# add some random valences to each edge
# valence can be 1 (positive tie) or -1 (negative tie)
# the 'valence' edge attribute is required
#E(g)$valence <- sample(c(-1, 1), length(E(g)), replace=T)
V(g)$name <- letters[1:N]

pii(g)


# there are 18 nodes
x <- read.csv('sampson.txt',stringsAsFactors=F)
N <- 18
networkNames <- c("samplk1","samplk2","samplk3","sampdlk","sampes","sampdes","sampin","sampnin","samppr","sampnpr")
x$Node <- trim(x$Node)

samp <- list()
for(i in 1:length(networkNames)) {
  y <- x[(1:18) + ((i - 1) * 18),]
  rownames(y) <- y$Node
  y <- as.matrix(y[,-1])
  samp[[i]] <- y
}
names(samp) <- networkNames


g <- igraph::graph.adjacency(samp[['samplk3']])
plotty(g)

# indegree for everything...
do.call('rbind', lapply(samp, function(x) {
  g <- graph.adjacency(x)
  degree(g, mode="in")
}))

# The political networks paper uses Esteem and Disesteem (sampes and sampdes)

g.pos <- graph.adjacency(samp[['sampes']])
g.neg <- graph.adjacency(samp[['sampdes']])

samp.el <- rbind(data.table(melt(samp[['sampes']]))[value != 0,][,valence := 1],
                 data.table(melt(samp[['sampdes']]))[value != 0,][,valence :=-1])
setnames(samp.el, c('Var1','Var2'),c('source','target'))



g <- graph.edgelist(as.matrix(samp.el[,c('source','target'),with=F]), directed = F)
E(g)$valence <- samp.el$valence
E(g)$weight <- samp.el$value
g <- simplify(g, edge.attr.comb = 'min')
g.samp <- g

g.pos <- graph.edgelist(as.matrix(samp.el[valence > 0,c('source','target'),with=F]), directed = F)
g.neg <- graph.edgelist(as.matrix(samp.el[valence < 0,c('source','target'),with=F]), directed = F)
E(graph.intersection(g.pos, g.neg))


### Testing PII ###############################################################

### The small example from the paper
# g <- as.undirected(graph.empty()) + vertices(letters[1:9])
# g <- g + path('i','e','a','c','g')
# g <- g + path('h','d','a','b','f')
# E(g)$valence <- 1
# E(g)[vid('a') %--% vid('c')]$valence <- -1
# E(g)[vid('a') %--% vid('b')]$valence <- -1

### An example (figure 2) from the paper
getThinkGraph <- function() {
  el <- matrix(c('a','b','a','c','b','c', 'a','d', 'd','e', 'd','f', 'd','g', 'g','e'), ncol=2, byrow=T)
  g <- graph.edgelist(el, directed = F)
  E(g)$valence <- 1
  E(g)[1 %--% 2]$valence <- -1
  E(g)[1 %--% 3]$valence <- -1
  E(g)[5 %--% 7]$valence <- -1
  E(g)$color <- ifelse(E(g)$valence == -1, "red", "black")
  g
}

getThinkGraph2 <- function() {
  el <- matrix(c('a','b', 'a','c', 'b','c', 'b','d', 'c','d', 'd','e', 'e','f'), ncol=2, byrow=T)
  g <- graph.edgelist(el, directed = F)
  E(g)$valence <- 1
  E(g)[5 %--% 6]$valence <- -1
  E(g)$color <- ifelse(E(g)$valence == -1, "red", "black")
  g
}
g <- getThinkGraph2()
plot(g)
pii(g, triadic=T)

cliques(g, min = 3, max = 3)
triad_holder = cliques(g, min = 3, max = 3)
node_distance_holder <- shortest.paths(g, V(g))

triadCalcs <- function(g){
  triad_table <- data.table(triadID = character(), nodeID = character(), direction = character(), valence = integer(),
                          distance = numeric(), nodeNum = integer())
  for(i in 1:length(triad_holder)){
    node1 = V(g)[triad_holder[[i]][1]]$name
    node2 = V(g)[triad_holder[[i]][2]]$name
    node3 = V(g)[triad_holder[[i]][3]]$name
    triname = paste(node1, node2, node3, sep = "-")
    edge1_2 = paste(node1, node2, sep = "-")
    edge1_3 = paste(node1, node3, sep = "-")
    edge2_3 = paste(node2, node3, sep = "-")
    for(number in 1:length(V(g))){
      refnode = V(g)[number]$name
      dis1 = edge.distance(g)[edge1_2, refnode]
      dis2 = edge.distance(g)[edge1_3, refnode]
      dis3 = edge.distance(g)[edge2_3, refnode]
      if(dis1 == dis2 & dis1 == dis3){
        nodedis1_2 = node_distance_holder[node1, node2]
        nodedis1_3 = node_distance_holder[node1, node3]
        nodedis2_3 = node_distance_holder[node2, node3]
        if(nodedis1_2 !=  nodedis1_3){
          if(nodedis1_2 != nodedis2_3){
            ref = edge1_2
          }else{ref = edge1_3}
        }else{ref = edge2_3}
        if(ref==edge1_2){truev = valence_set[1]}
        if(ref==edge1_3){truev = valence_set[3]}
        if(ref==edge2_3){truev = valence_set[2]}
        row <- data.table(triadID = triname, nodeID = refnode, direction = 'IN ', valence = truev,
                          distance = edge.distance(g)[edge1_2, refnode], nodeNum = as.integer(V(g)[refnode]))
        triad_table <- rbind(triad_table, row, use.names=T)
      }

      else{
        if(dis1 != dis2 ){
          if(dis1 != dis3){
            diff = edge1_2
          }else{diff = edge1_3}
        }else{diff = edge2_3}
        p <- c(triad_holder[[i]], triad_holder[[i]][1])
        valence_set <- E(g, path = p)$valence
        #print(valence_set)
        if(diff==edge1_2){truev = valence_set[1]}
        if(diff==edge1_3){truev = valence_set[3]}
        if(diff==edge2_3){truev = valence_set[2]}
        row <- data.table(triadID = triname, nodeID = refnode, direction = 'OUT', valence = truev,
                          distance = edge.distance(g)[diff, refnode],  nodeNum = as.integer(V(g)[refnode]))
        triad_table <- rbind(triad_table, row, use.names=T)
      }
    }
  }
  return(triad_table)
}
triad_table <- triadCalcs(g)

(regularPII <- function(pii.beta = -0.8){
  if(is.null(e.dist)) {
    e.dist <- edge.distance(g)
  }
  e.dist <- matrix(as.integer(e.dist), nrow=nrow(e.dist)) # convert to an integer matrix
  max.distance <- max(e.dist)
  max.degree <- max(degree(g, mode='total'))
  valence <- E(g)$valence
  pii.x <- (log(2) - log(abs(pii.beta))) / log(max.degree)
  #numEdge <- length(E(g))
  #numNode <- length(V(g))
  print(numEdge)
  print(numNode)

  #nNegCount = c()
  #outNegCount = c()
  #inPosCount = c()
  #outPosCount = c()

  piiBetaVector = c()
  piIndex = numeric(numNode)
  for(a in 0:max.distance){
    piiBetaVector <- c(piiBetaVector, pii.beta^a)
  }
  for(i in 1:numNode){
    negEdgeCount = numeric(max.distance+1)
    posEdgeCount = numeric(max.distance+1)
    for(j in 1:numEdge){
      if(valence[j] < 0){
        negEdgeCount[edge.distance(g)[j,i]+1] <- negEdgeCount[edge.distance(g)[j,i]+1]+1
      }else{
        posEdgeCount[edge.distance(g)[j,i]+1] <- posEdgeCount[edge.distance(g)[j,i]+1]+1
      }
    }
    for(k in 1:(max.distance+1)){
      piIndex[i] = piIndex[i] + (piiBetaVector[k] * (posEdgeCount[k]^pii.x - negEdgeCount[k]^pii.x))
    }
    #cat('Node',i, ' counts \n')
#     print(negEdgeCount)
#     print(posEdgeCount)
     print(piIndex)
  }
})()

(triadicPII <- function(g, pii.beta = -0.8, pii.delta = 0.1, e.dist = NULL){
  if(is.null(e.dist)) {
    e.dist <- edge.distance(g)
  }
  e.dist <- matrix(as.integer(e.dist), nrow=nrow(e.dist)) # convert to an integer matrix
  max.distance <- max(e.dist)
  max.degree <- max(degree(g, mode='total'))
  valence <- E(g)$valence
  pii.x <- (log(2) - log(abs(pii.beta))) / log(max.degree)
  numEdge <- length(E(g))
  numNode <- length(V(g))
#   print(numEdge)
#   print(numNode)

  piiBetaVector = c()
  piIndex = numeric(numNode)
  for(a in 0:max.distance){
    piiBetaVector <- c(piiBetaVector, pii.beta^a)
  }
  for(i in 1:numNode){
    negEdgeCount = numeric(max.distance+1)
    posEdgeCount = numeric(max.distance+1)
    inNegCount = numeric(max.distance+1)
    outNegCount = numeric(max.distance+1)
    inPosCount = numeric(max.distance+1)
    outPosCount = numeric(max.distance+1)
    for(j in 1:numEdge){
      if(valence[j] < 0){
        negEdgeCount[edge.distance(g)[j,i]+1] <- negEdgeCount[edge.distance(g)[j,i]+1]+1
      }else{
        posEdgeCount[edge.distance(g)[j,i]+1] <- posEdgeCount[edge.distance(g)[j,i]+1]+1
      }
    }
    for(r in 1:length(triad_table$valence)){
      currentValence = triad_table$direction[r]
      if(triad_table$nodeID[r] == V(g)$name[i]){
        if(triad_table$direction[r] == "OUT"){
          if(triad_table$valence[r] < 0){
            outNegCount[triad_table$distance[r]+1] <- outNegCount[triad_table$distance[r]+1]+1
          }else{
            outPosCount[triad_table$distance[r]+1] <- outPosCount[triad_table$distance[r]+1]+1
          }
        }else{
          if(triad_table$valence[r] < 0){
            inNegCount[triad_table$distance[r]+1] <- inNegCount[triad_table$distance[r]+1]+1
          }else{
            inPosCount[triad_table$distance[r]+1] <- inPosCount[triad_table$distance[r]+1]+1
          }
        }
      }
    }
    for(k in 1:(max.distance+1)){
      piIndex[i] = piIndex[i] + (piiBetaVector[k] * (posEdgeCount[k]^pii.x - negEdgeCount[k]^pii.x + pii.delta*(
        outNegCount[k]^pii.x - outPosCount[k]^pii.x + inNegCount[k]^pii.x - inPosCount[k]^pii.x)))
    }
#     cat('Node',i, ' counts \n')
#     print(inNegCount)
#     print(outNegCount)
#     print(outPosCount)
#     print(inPosCount)
#     #print(negEdgeCount)
#     #print(posEdgeCount)
#     print(piIndex)
  }
  return(piIndex)
})(g)

piis <- data.table(node = character(), pii = numeric(), delta = numeric())
for(d in seq(0,1,by=0.1)) {
  cat('Delta: ', d, '\n')
  td <- triadicPII(g, pii.delta = d)
  piis <- rbind(piis, data.table(node = V(g)$name, pii = td, delta = d))
}

ggplot(piis) + geom_line(aes(x=delta, y=pii, group=node, color=node))


triadicPII(g, pii.delta=0.1)

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
plot(randomGraph())


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


piis <- data.table(nd=character(), pii=numeric(), beta=numeric())
for(b in seq(-1, 0, by=0.001)) {
  piis <- rbind(piis, data.table(nd = V(g)$name, pii=pii(g, pii.beta=b), beta = b))
}
ggplot(piis) + geom_line()


#centralization.degree(g)$centralization
#vcount(g)
#ecount(g)
#graph.density(g)
#length(E(g)[valence == -1])
#length(E(g)[valence == 1])
#average.path.length(g)
###var(pii)
#beta.sequence <- seq(-1, -0.1, by=0.1)
#all.pii <- data.table(node=character(), pii.value = numeric(), beta = numeric())
#g <- all.ring.graphs[[1]]
#g.ed <- edge.distance(g)
#for(b in beta.sequence) {
#  g.pii <- pii(g,e.dist = g.ed, pii.beta = b)
#  td <- data.table(node=1:vcount(g), pii.value=as.numeric(g.pii), degree = degree(g), beta=b)
#  all.pii <- rbind(all.pii, td)
#}

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
#md(all.ring.graphs[[1]])
gp <- do.call('rbind', mclapply(all.graphs, graphData, mc.cores=5))
ggplot(gp, aes(y=rankCor, x=diameter)) + geom_point() + geom_smooth(method='lm')

### Random graphs - Watts Strogatz ############################################
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
randomGraph <- function(i=0){
  N <- sample(10:100, 1)
  graph <- watts.strogatz.game(1, N, 3, 0.05)
  n <- sample(1:80, 1) / 100
  E(graph)$valence <- sample(c(-1, 1), ecount(graph), replace = T, prob = c(n, 1-n))
  E(graph)$color <- ifelse(E(graph)$valence == -1, "red", "black")
  V(graph)$name <- letterNames(vcount(graph))
  graphid <- paste0('watts-strogatz ', '#',i, sep='')
  graph <- set.graph.attribute(graph, 'graphid', graphid)
  graph
}
igraph.options(vertex.color = '#FFFFFFAA', vertex.alpha = 0.5, vertex.size = 15, vertex.frame.color = 'grey60',
               vertex.label.family = 'Ubuntu', vertex.label.cex = 0.7, vertex.label.color = 'blue')
g <- randomGraph()
plot(g)

rand.graphs = list()
for(i in 1:10000){
  rand.graphs[[i]] <- randomGraph(i)
}

gp <- do.call('rbind', mclapply(rand.graphs, graphData, mc.cores=5))
ggplot(gp, aes(y=rankCor, x=avgPathLength)) + geom_point() + geom_smooth(method='lm')
ggplot(gp, aes(y=rankCor, x=propNegEdge)) + geom_point() + geom_smooth(method='lm')
summary(lm1 <- lm(rankCor ~ meanTrans + avgPathLength + degCentralization + propNegEdge + density + modularity, data=gp))
summary(lm2 <- lm(rankCor ~ avgPathLength + propNegEdge + avgPathLength:numNegEdge, data=gp))
summary(lm2 <- lm(rankCor ~ modularity, data=gp))

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


library(parallel)
all.graphs <- list()
for(i in 3:10) {
  all.graphs <- c(all.graphs, getAllRingGraphs(i))
  all.graphs <- c(all.graphs, getAllRingGraphs(i, chain=T))
}
nd <- do.call('rbind', mclapply(all.graphs, nodeData, mc.cores=5))

# cnd <- do.call('rbind', lapply(getAllRingGraphs(3, chain = T), nodeData))
# for(i in 4:10){
#   cnd <- rbind(cnd, do.call('rbind', lapply(getAllRingGraphs(i, chain = T), nodeData)))
# }


#' table 3
#' graph level
#'


nd3 <- group_by(nd, graphid, beta) %>%
  summarize(meanPII = mean(pii.value), varPII = var(pii.value))


all.graph.name <- unique(nd$graphid)
k <- 31
pdf(paste0('catalog.pdf'))
for(k in 1:length(all.graph.name))  {
  plot(all.graphs[[all.graph.name[k]]])
  p <- ggplot(nd[graphid == all.graph.name[k]]) +
    geom_line(size=3, alpha=0.6, position=position_jitter(height=0.02, width=0.02),
              aes(x=beta, y=pii.value, group=nodeid, color=nodeid))
  print(p)
}
dev.off()


#centralization.degree(g)$centralization
#vcount(g)
#ecount(g)
#graph.density(g)
#length(E(g)[valence == -1])
#length(E(g)[valence == 1])
#average.path.length(g)
###var(pii)
beta.sequence <- seq(-1, -0.1, by=0.1)
delta.sequence <- seq(0, 1, by=0.1)
all.pii <- data.table(node=character(), pii.value = numeric(), beta = numeric(), delta = numeric())
g.ed <- edge.distance(g)
for(i in 1:length(beta.sequence)) {
  for(j in 1:length(delta.sequence)) {
    cat('beta: ', i, ' - delta: ', j,'\r')
    samp.pii <- pii(g,e.dist = g.ed, pii.beta = beta.sequence[i], triadic = T, pii.delta = delta.sequence[j])
    all.pii <- rbind(all.pii,
                     data.table(node=names(samp.pii), pii.value=as.numeric(samp.pii),
                                beta=beta.sequence[i], delta=delta.sequence[j]))
  }
}







all.pii[, is.neg := ifelse(pii.value < 0, "neg", "pos")]
all.pii[, pii.value := rescale(pii.value)]
ggplot(all.pii, aes(x=beta, y=delta, fill=pii.value)) + geom_tile() + facet_wrap(~ node) +
  scale_x_continuous('Beta') + scale_y_continuous('Delta') +
  scale_fill_gradient2()

all.pii.agg <- all.pii[,list(m.pii = mean(pii.value, na.rm=T)), by=list(beta, delta)]
ggplot(all.pii.agg, aes(x=beta, y=delta, fill=m.pii)) + geom_tile() +
  scale_x_continuous('Beta') + scale_y_continuous('Delta') +
  scale_fill_gradient2(low = "white", mid="brown", high = "steelblue")

ggplot(all.pii.agg) + geom_point(size=15, aes(x = beta, y=delta, color=m.pii))


### Compute an auto-correlation matrix

d <- 0.5
all.pii.d <- all.pii[delta == d]
psi.corr <- data.table(beta1 = numeric(), beta2 = numeric(), psi = numeric())
for(b1 in unique(all.pii.d$beta)) {
  for(b2 in unique(all.pii.d$beta)) {
    psi <- cor(all.pii.d[beta == b1]$pii.value, all.pii.d[beta == b2]$pii.value)
    psi.corr <- rbind(psi.corr, data.table(beta1=b1, beta2=b2, psi=psi))
  }
}

ggplot(psi.corr, aes(x=beta1, y=beta2, fill=psi)) + geom_tile()

psi.mat <- matrix(psi.corr[order(beta1, beta2)]$psi, nrow=10)

eigen(psi.mat)


### benchmark different versions
e.dist <- edge.distance(g)
benchmark(pii(g), pii2(g), pii(g, e.dist = e.dist), pii2(g, e.dist = e.dist), pii2(g, e.dist = e.dist, triadic = T), replications=10, order='relative')
benchmark(edge.distance(g), edge.distance2(g), replications = 10, order = 'relative')


g.lo <- layout.kamada.kawai(g)

x <- pii2(g, pii.beta = -0.8)
plot(g, layout=g.lo, vertex.size=rescale(x)*50, vertex.color='SkyBlue2',
     vertex.frame.color='SkyBlue2', edge.arrow.size=0.5,
     vertex.label.family='Open Sans', vertex.label.cex=2)



