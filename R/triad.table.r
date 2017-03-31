#' Triad Table
#'
#' Creates a table that contains triadID - the three nodes that make up the triad, nodeID - the character name of the node,
#' direction of the triad from the node, the valence of the triad from the node's perspective, distance of the triad from the node,
#' and nodeNum - the numeric value of the node
#'
#' @param g An igraph graph
#' @export
triad.table <- function(g) {
  if(is.null(V(g)$name)){V(g)$name <- 1:vcount(g)}
  #e.dist <- edge.distance(g)
  node_distance_holder <- shortest.paths(g, V(g))
  #e.dist.l <- list(ncol(e.dist))
  #nms <- strsplit(rownames(e.dist), "-")
  # #for(j in 1:ncol(e.dist)) {
  #   m1 <- matrix(NA, nrow = vcount(g), ncol = vcount(g))
  #   rownames(m1) <- colnames(e.dist)
  #   colnames(m1) <- colnames(e.dist)
  #
  #   for(i in 1:nrow(e.dist)) {
  #     x <- nms[[i]]
  #     m1[x[1], x[2]] <- e.dist[i, j]
  #     m1[x[2], x[1]] <- e.dist[i, j]
  #   }
  #   e.dist.l[[j]] <- m1
  # }
  e.dist.l <- edge.distance(g, lookup.mat = T)
  triads <- cliques(g, min = 3, max= 3)
  triads <- do.call('rbind', lapply(triads, function(x) { as.integer(x) }))
  dimNum <- vcount(g)
  if(is.null(triads)){stop("There are no triads in the graph")}
  triad_table <- triadTable(edgeDistance = e.dist.l, shortPaths = node_distance_holder, triads=triads, vertices=V(g), edgevalence = E(g)$valence)
  return(triad_table)
}
