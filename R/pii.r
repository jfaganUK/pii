#' Political Independence Index
#'
#' Returns a numeric vector of PII scores for each vertex in network.
#'
#' @useDynLib pii
#' @importFrom Rcpp sourceCpp
#' @param g An igraph graph
#' @param pii.beta Should the vertex and edge names be added to the rows and columns of the matrix.
#' @param e.dist (optional) an edge.distance matrix, if calculated ahead of time.
#' @export
#' @examples
#' pii(g)

pii <- function(g, pii.beta = -0.8, e.dist = NULL, triadic = F, pii.delta = 0.1) {
  if(!("igraph" %in% class(g))) {
    stop("The graph object must be an igraph object.")
  }
  # TODO: Instead of calculating per component. Give the users an option to calculate for the whole network
  # Also, if they don't want to calculate for the "whole-network" then we still want PII values for
  # each node, but just run pii for each component. If you are doing the whole network though, calculate
  # 'x' first and use the same 'x' value for each run.
  if((cl <- igraph::clusters(g))$no > 1) {
    warning("Disconnected graph. Using only the first giant component.")
    g <- induced.subgraph(g, which(cl$membership == which.max(cl$csize)))
  }
  if(!("valence" %in% names(edge.attributes(g)))) {
    warning("Valence attribute not found. Assuming all ties are positive.")
    E(g)$valence <- 1
  }
  if(length(E(g)) < 1) {
    stop("Only one edge in the network.")
  }
  if(is.null(e.dist)) {
    e.dist <- edge.distance(g)
  }
  e.dist <- matrix(as.integer(e.dist), nrow=nrow(e.dist)) # convert to an integer matrix
  max.distance <- max(e.dist)
  max.degree <- max(degree(g, mode='total'))
  edgevalence <- E(g)$valence
  pii.x <- (log(2) - log(abs(pii.beta))) / log(max.degree)
  if(triadic) {
    triad_table <- triadCalcs(g)
    x <- piiTriadicCalc(e.dist, edgevalence, pii.beta, pii.x, max.distance, triad_table, pii.delta)
  } else {
    x <- piiCalc(e.dist, edgevalence, pii.beta, pii.x, max.distance)
  }
  names(x) <- V(g)$name
  x
}
