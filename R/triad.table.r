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
  if(is.null(triads)){stop("There are no triads in the graph")}

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

  b$nodeDistA <- n.dist[as.matrix(b[, c('triadNodeA', 'focalNode')])]
  b$nodeDistB <- n.dist[as.matrix(b[, c('triadNodeB', 'focalNode')])]
  b$nodeDistC <- n.dist[as.matrix(b[, c('triadNodeC', 'focalNode')])]

  # direction, 1 is inward, 0 is outward
  b <- mutate(b, direction = ifelse(distEdgeAB == distEdgeAC & distEdgeAB == distEdgeBC,
                                    ifelse(nodeDistA == nodeDistB & nodeDistA == nodeDistC, 2, 1), 0))
  b$closingNode1 <- NA
  b$closingNode2 <- NA

  #lookup matrix for valence of edges between two lookup nodes
  v.mat <- matrix(NA, ncol = vcount(g), nrow = vcount(g))
  v.mat[get.edgelist(g, names = F)] <- E(g)$valence

  # with outward facing triads, the closing edge is always the furthest away
  farthestEdge <- apply(b[b$direction == 0,c('distEdgeAB', 'distEdgeAC', 'distEdgeBC')], 1, which.max)

  m <- matrix(c(1, 2, 1, 3, 2, 3), byrow = T, ncol = 2)
  m <- m[farthestEdge, ]
  m <- cbind(1:nrow(m), m)

  b1 <- b[b$direction == 0, c('triadNodeA', 'triadNodeB', "triadNodeC")]
  b[b$direction == 0, "closingNode1"] <- b1[m[, c(1, 2)]]
  b[b$direction == 0, "closingNode2"] <- b1[m[, c(1, 3)]]

  # with inward facing triads, the closing edge is the one between the two closest nodes
  b1 <- b[b$direction == 1, c('triadNodeA', 'triadNodeB', "triadNodeC")]
  x <- apply(b[b$direction == 1, c("nodeDistA", "nodeDistB", "nodeDistC")], 1, function(x) { which(x == min(x)) }) %>%
     t %>% as.matrix %>% cbind(1:nrow(.), .)
  b[b$direction == 1, 'closingNode1'] <- b1[x[, c(1, 2)]]
  b[b$direction == 1, 'closingNode2'] <- b1[x[, c(1, 3)]]

  # ambiguous triads are replicated and given the valence of each edge
  b1 <- b[b$direction == 2, ]
  b <- b[b$direction != 2,]

  b1a <- b1
  b1a$closingNode1 <- b1a$triadNodeA
  b1a$closingNode2 <- b1a$triadNodeB

  b1b <- b1
  b1b$closingNode1 <- b1b$triadNodeA
  b1b$closingNode2 <- b1b$triadNodeC

  b1c <- b1
  b1c$closingNode1 <- b1c$triadNodeB
  b1c$closingNode2 <- b1c$triadNodeC

  b <- rbind(b, b1a, b1b, b1c)

  # add valence
  b$closeEdgeValence <- v.mat[b[,c("closingNode1", "closingNode2")] %>% as.matrix]
  b$closeEdgeDist <- e.dist.l[b[,c("closingNode1", "closingNode2", "focalNode")] %>% as.matrix]

  return(b)
}
