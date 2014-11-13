#' Triadic Edges
#'
#' Identifies which edges in a network are part of a triad.
#'
#' @param g An igraph graph
#' @export
#' @examples
#' get.edge.network(g)

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
