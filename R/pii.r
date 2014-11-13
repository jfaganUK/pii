#' Political Independence Index
#'
#' Returns a numeric vector of PII scores for each vertex in the network.
#'
#' @param g An igraph graph
#' @param pii.beta Should the vertex and edge names be added to the rows and columns of the matrix.
#' @param e.dist (optional) an edge.distance matrix, if calculated ahead of time.
#' @export
#' @examples
#' pii(g)

pii <- function(g, pii.beta = -0.8, e.dist = NULL) {
  if(class(g) != 'igraph') {
    stop('g should be an igraph object.')
  }
  if(is.null(e.dist)) {
    e.dist <- edge.distance(g)
  }

  e.dist <- data.table(melt(e.dist))
  setnames(e.dist, c("Var1","Var2", "value"), c('enm','node', 'distance'))
  e.dist$enm <- as.character(e.dist$enm)
  e.dist$node <- as.character(e.dist$node)

  max.distance <- max(e.dist$distance, na.rm=T)
  max.degree <- max(degree(g, mode='total'))

  # calulate the beta ^ i ahead of time, I think, maybe, this saves time, I don't know.
  pii.beta.exp <- numeric(max.distance+1)
  for(i in 0:max.distance) {
    pii.beta.exp[i+1] <- pii.beta ^ (i)
  }

  pii.x <- (log(2) - log(abs(pii.beta))) / log(max.degree)

  pos.ties <- matrix(0, nrow=vcount(g), ncol=max.distance+1)
  rownames(pos.ties) <- V(g)$name
  colnames(pos.ties) <- 0:max.distance

  neg.ties <- matrix(0, nrow=vcount(g), ncol=max.distance+1)
  rownames(neg.ties) <- V(g)$name
  colnames(neg.ties) <- 0:max.distance

  g.el <- data.table(get.edgelist(g))
  setnames(g.el,c("V1","V2"), c("source", "target"))
  g.el$enm <- paste(g.el[,source],g.el[,target],sep="-")
  g.el$source <- as.character(g.el$source)
  g.el$target <- as.character(g.el$target)
  g.el$valence <- E(g)$valence

  g.el.dist <- merge(e.dist, g.el, by='enm', all.x=T)
  count.distance <- g.el.dist[, .N, by=list(node, valence, distance)]

  pi.index <- numeric(vcount(g))
  nms <- V(g)$name
  names(pi.index) <- nms
  for(i in 1:length(pi.index)) {
    count.distance1 <- count.distance[node == nms[i]]
    for(j in 0:max.distance) {
      Pn <- count.distance1[distance == j & valence == 1, N] ^ pii.x
      Nn <- count.distance1[distance == j & valence == -1, N] ^ pii.x
      Pn <- ifelse(length(Pn) == 0, 0, Pn)
      Nn <- ifelse(length(Nn) == 0, 0, Nn)
    }
  }
  pi.index
}
