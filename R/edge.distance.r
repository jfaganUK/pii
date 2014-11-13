#' Edge Distance
#'
#' Calculates the distance from each node in the network to each edge in the network.
#'
#' @param g An igraph graph
#' @export
#' @examples
#' get.edge.network(g)

edge.distance <- function(g, edge.network = NULL, add.names=T) {
  if(is.null(edge.network)) {
    edge.network <- get.edge.network(g)
  }
  ns <- V(g)$name
  max.node.id <- vcount(g)
  enm <- (max.node.id + 1) : (max.node.id + ecount(g))
  g.dist <- shortest.paths(graph.edgelist(edge.network))
  g.dist <- (g.dist[enm, 1:max.node.id] - 1) / 2
  if(add.names) {
    if("matrix" %in% class(g.dist)) {
      colnames(g.dist) <- ns
    } else {
      names(g.dist) <- ns
    }
    g.el <- get.edgelist(g)
    rownames(g.dist) <- paste(g.el[,1],g.el[,2],sep="-")
  }
  g.dist
}
