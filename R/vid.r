#' Vertex ID
#'
#' Returns the vertex ID when given a name and a graph.
#'
#' @param g An igraph graph
#' @param x A vertex name (character)
#' @export
#' @examples
#' vid(g, x="Billy")

vid <- function(g, x) which(V(g)$name == x)
