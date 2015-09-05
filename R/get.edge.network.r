#' Get Edge Network
#'
#' Returns a modfied network where each edge in the original network is split in two and a new node is inserted in between and given the label of the original edge.
#'
#' @param g An igraph graph
#' @export
#' @examples
#' g <- graph.ring(10)
#' get.edge.network(g)

get.edge.network <- function(g) {
  require(igraph)
  el <- get.edgelist(g, names = F)
  el <- matrix(as.integer(el), nrow=nrow(el)) # convert to an integer matrix
  oel <- getEdgeNetworkCalc(el)
  return(oel)
}
