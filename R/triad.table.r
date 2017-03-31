#' Triad Table
#'
#' Creates a table that contains triadID - the three nodes that make up the triad, nodeID - the character name of the node,
#' direction of the triad from the node, the valence of the triad from the node's perspective, distance of the triad from the node,
#' and nodeNum - the numeric value of the node
#'
#' @param g An igraph graph
#' @useDynLib pii
#' @export
triad.table <- function(g) {
  if(is.null(V(g)$name)){V(g)$name <- 1:vcount(g)}
  n.dist <- shortest.paths(g, V(g))
  e.valence <- cbind(get.edgelist(g, names = F), E(g)$valence)

  # 3d edge distance matrix of edge, edge, focal.node
  e.dist.l <- edge.distance(g, lookup.mat = T)

  # enumerate all the triads
  triads <- cliques(g, min = 3, max= 3)
  triads <- do.call('rbind', lapply(triads, function(x) { as.integer(x) }))

  # start the triad table
  # create a copy of the triads for each node, stack it all up
  b <- do.call(rbind, replicate(vcount(g), triads, simplify=F)) %>%
    cbind(rep(1:vcount(g), each = nrow(triads))) %>%
    data.frame %>%
    rename(triadNodeA = X1, triadNodeB = X2, triadNodeC = X3, focalNode = X4)

  # use the two nodes in each edge and the focal node as a lookup index in the e.dist.l table
  b$distEdgeAB <- e.dist.l[as.matrix(b[,c('triadNodeA','triadNodeB', 'focalNode')])]
  b$distEdgeAC <- e.dist.l[as.matrix(b[,c('triadNodeA','triadNodeC', 'focalNode')])]
  b$distEdgeBC <- e.dist.l[as.matrix(b[,c('triadNodeB','triadNodeC', 'focalNode')])]

  # direction, 1 is inward, 0 is outward
  b <- mutate(b, direction = ifelse(distEdgeAB == distEdgeAC & distEdgeAB == distEdgeBC, 1, 0))
  b$closingEdge1 <- NA
  b$closingEdge2 <- NA

  # with outward facing triads, the closing edge is always the furthest away
  farthestEdge <- apply(b[b$direction == 0,c('distEdgeAB', 'distEdgeAC', 'distEdgeBC')], 1, which.max)
  farthestEdge

  # with inward facing triads, the closing edge is the one between the two closest nodes
  b$nodeDistA <- n.dist[as.matrix(b[, c('triadNodeA', 'focalNode')])]
  b$nodeDistB <- n.dist[as.matrix(b[, c('triadNodeB', 'focalNode')])]
  b$nodeDistC <- n.dist[as.matrix(b[, c('triadNodeC', 'focalNode')])]





  dimNum <- vcount(g)
  if(is.null(triads)){stop("There are no triads in the graph")}
  triad_table <- triadTable(edgeDistance = e.dist.l, shortPaths = n.dist, triads=triads, vertices=V(g), edgevalence = E(g)$valence)
  return(triad_table)
}
