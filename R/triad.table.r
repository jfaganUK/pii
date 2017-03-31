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
  n.dist <- shortest.paths(g, V(g))
  e.dist.l <- edge.distance(g, lookup.mat = T)
  triads <- cliques(g, min = 3, max= 3)
  triads <- do.call('rbind', lapply(triads, function(x) { as.integer(x) }))
  dimNum <- vcount(g)
  if(is.null(triads)){stop("There are no triads in the graph")}
  triad_table <- triadTable(edgeDistance = e.dist.l, shortPaths = n.dist, triads=triads, vertices=V(g), edgevalence = E(g)$valence)
  return(triad_table)
}
