# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

getEdgeNetworkCalc <- function(inputEdgeList) {
    .Call('pii_getEdgeNetworkCalc', PACKAGE = 'pii', inputEdgeList)
}

#' Calculate the PII scores after the edge distances are calculated
#' 
#' @param edgeDistance An integer matrix of distances from the nodes to the edges
#' @param valence The postive / negative valence of each of the edges.
#' @param piiBeta The beta attenuation value (ideally -1.0 to -0.01)
#' @param piiX The normalization factor based on the maximum degree.
#' @param maxDistance The maximum distance in the edgeDistance matrix (or some user-provided max distance for very large networks)
#' @export
piiCalc <- function(edgeDistance, valence, piiBeta, piiX, maxDistance) {
    .Call('pii_piiCalc', PACKAGE = 'pii', edgeDistance, valence, piiBeta, piiX, maxDistance)
}

piiTriadicCalc <- function(edgeDistance, valence, piiBeta, piiX, maxDistance, edgeTriadic, piiDelta) {
    .Call('pii_piiTriadicCalc', PACKAGE = 'pii', edgeDistance, valence, piiBeta, piiX, maxDistance, edgeTriadic, piiDelta)
}
