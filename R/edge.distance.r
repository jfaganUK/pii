#' Edge Distance
#'
#' Returns a rectangular matrix where the rows are the edges of the graph, and the columns are the vertices of the graph.
#' The cells of the matrix is the distance from the vertex to the edge. The edges directly incident on the vertex are distance 0.
#' It is used in the calculation of PII.
#'
#' @param g An igraph graph
#' @param add.names Should the vertex and edge names be added to the rows and columns of the matrix.
#' @export
#' @examples
#' edge.distance(g)

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
