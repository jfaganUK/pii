
vid <- function(g, x) which(V(g)$name == x)

edge.distance <- function(g, add.names=T) {
  g.el <- get.edgelist(g, names=F)
  ns <- V(g)$name
  
  max.node.id <- max(g.el)
  enm <- (max.node.id + 1) : (max.node.id + nrow(g.el))
  n <- nrow(g.el) * 2
  new.el <- matrix(0, nrow=n, ncol=2)
  
  for(i in 1:nrow(g.el)) {
    # from v1 to the edge
    new.el[(i*2)-1, 1] <- g.el[i, 1]
    new.el[(i*2)-1, 2] <- enm[i]
    
    # from the edge to v2
    new.el[(i*2),   1] <- enm[i]
    new.el[(i*2),   2] <- g.el[i, 2]
  }
  
  g.dist <- shortest.paths(graph.edgelist(new.el))
  g.dist <- (g.dist[enm, 1:max.node.id] - 1) / 2
  if(add.names) {
    colnames(g.dist) <- ns
    g.el <- get.edgelist(g)
    rownames(g.dist) <- paste(g.el[,1],g.el[,2],sep="-")
  }
  g.dist
}

pii <- function(g, pii.beta = -0.8, e.dist = NULL) {
  if(!require(igraph)) {
    stop('Requires igraph package')
  }
  if(!require(data.table)) {
    stop('Requires the data.table package')
  }
  if(!require(reshape2)) {
    stop('Requires the reshape2 package')
  }
  if(class(g) != 'igraph') {
    stop('g should be an igraph object.')
  }
  if(is.null(e.dist)) {
    e.dist <- edge.distance(g)
  }
  
  e.dist <- data.table(melt(e.dist))
  setnames(e.dist, c("Var1","Var2", "value"), c('enm','node', 'distance'))
  e.dist$enm <- as.character(e.dist$enm)
  e.dist$node <- as.character(e.dist$node)
  
  max.distance <- max(e.dist$distance, na.rm=T)
  max.degree <- max(degree(g, mode='total'))
  
  # calulate the beta ^ i ahead of time, I think, maybe, this saves time, I don't know.
  pii.beta.exp <- numeric(max.distance+1)
  for(i in 0:max.distance) {
    pii.beta.exp[i+1] <- pii.beta ^ (i)
  }
  
  pii.x <- (log(2) - log(abs(pii.beta))) / log(max.degree)
  
  pos.ties <- matrix(0, nrow=vcount(g), ncol=max.distance+1)
  rownames(pos.ties) <- V(g)$name
  colnames(pos.ties) <- 0:max.distance
  
  neg.ties <- matrix(0, nrow=vcount(g), ncol=max.distance+1)
  rownames(neg.ties) <- V(g)$name
  colnames(neg.ties) <- 0:max.distance
  
  g.el <- data.table(get.edgelist(g))
  setnames(g.el,c("V1","V2"), c("source", "target"))
  g.el$enm <- paste(g.el[,source],g.el[,target],sep="-")
  g.el$source <- as.character(g.el$source)
  g.el$target <- as.character(g.el$target)
  g.el$valence <- E(g)$valence
  
  g.el.dist <- merge(e.dist, g.el, by='enm', all.x=T)
  count.distance <- g.el.dist[, .N, by=list(node, valence, distance)]
  
  pi.index <- numeric(vcount(g))
  nms <- V(g)$name
  names(pi.index) <- nms
  for(i in 1:length(pi.index)) {
    count.distance1 <- count.distance[node == nms[i]]
    for(j in 0:max.distance) {
      Pn <- count.distance1[distance == j & valence == 1, N] ^ pii.x
      Nn <- count.distance1[distance == j & valence == -1, N] ^ pii.x
      Pn <- ifelse(length(Pn) == 0, 0, Pn)
      Nn <- ifelse(length(Nn) == 0, 0, Nn)
    }
  }
  pi.index
}

sourceCpp('pii.cpp')

pii2 <- function(g, pii.beta = -0.8, e.dist = NULL, triadic = F, pii.delta = 1, e.triads = NULL) {
  if(is.null(e.dist)) {
    e.dist <- edge.distance(g)    
  }
  e.dist <- matrix(as.integer(e.dist), nrow=nrow(e.dist)) # convert to an integer matrix
  max.distance <- max(e.dist)
  max.degree <- max(degree(g, mode='total'))
  valence <- E(g)$valence
  pii.x <- (log(2) - log(abs(pii.beta))) / log(max.degree)
  if(triadic) {
    if(is.null(e.triads)) {  
      e.triads <- triadic.edges(g)
    }
    x <- piiTriadicCalc(e.dist, valence, pii.beta, pii.x, max.distance, e.triads, pii.delta)
  } else {
    x <- piiCalc(e.dist, valence, pii.beta, pii.x, max.distance)    
  }
  names(x) <- V(g)$name
  x
}

sourceCpp('getEdgeNetwork.cpp')
get.edge.network <- function(g) {
  el <- get.edgelist(g, names = F)
  el <- matrix(as.integer(el), nrow=nrow(el)) # convert to an integer matrix
  oel <- getEdgeNetworkCalc(el)
  return(oel)
}

edge.distance2 <- function(g, edge.network = NULL, add.names=T) {
  if(is.null(edge.network)) {
    edge.network <- get.edge.network(g)
  }
  ns <- V(g)$name
  max.node.id <- vcount(g)
  enm <- (max.node.id + 1) : (max.node.id + ecount(g))
  g.dist <- shortest.paths(graph.edgelist(edge.network))
  g.dist <- (g.dist[enm, 1:max.node.id] - 1) / 2
  if(add.names) {
    colnames(g.dist) <- ns
    g.el <- get.edgelist(g)
    rownames(g.dist) <- paste(g.el[,1],g.el[,2],sep="-")
  }
  g.dist
}

triadic.edges <- function(g) {
  all.triads <- graph.get.subisomorphisms.vf2(g, graph.ring(3))
  convert.nl.to.el <- function(x) {
    matrix(c(x[1], x[2], x[1], x[3], x[2], x[3]), nrow=3, byrow=T)
  }
  triad.edges <- do.call('rbind', lapply(all.triads, 'convert.nl.to.el'))
  triad.edges <- unique(paste(triad.edges[,1],triad.edges[,2],sep="-"))
  g.el <- get.edgelist(g, names=F)
  g.el <- paste(g.el[,1],g.el[,2],sep="-")
  is.triadic <- g.el %in% triad.edges
  return(is.triadic)
}


