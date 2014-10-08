
sourceCpp('pii.cpp')

pii2 <- function(g, pii.beta = -0.8, e.dist = NULL) {
  if(is.null(e.dist)) {
    e.dist <- edge.distance(g)    
  }
  e.dist <- matrix(as.integer(e.dist), nrow=nrow(e.dist)) # convert to an integer matrix
  max.distance <- max(e.dist)
  max.degree <- max(degree(g, mode='total'))
  valence <- E(g)$valence
  pii.x <- (log(2) - log(abs(pii.beta))) / log(max.degree)
  x <- piiCalc(e.dist, valence, pii.beta, pii.x, max.distance)
  names(x) <- nds
  x
}
