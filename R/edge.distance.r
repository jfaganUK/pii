#' Edge Distance
#'
#' Calculates the distance from each node in the network to each edge in the network.
#'
#' @param g An igraph graph
#' @param edge.network A pre-calculated edge network from ```get.edge.network```
#' @export
#' @examples
#' get.edge.network(g)

edge.distance <- function(g, edge.network = NULL, lookup.mat = F, add.names=T) {
  if(is.null(edge.network)) {
    edge.network <- get.edge.network(g)
  }
  ns <- V(g)$name
  max.node.id <- vcount(g)
  enm <- (max.node.id + 1) : (max.node.id + ecount(g))
  g.dist <- shortest.paths(graph.edgelist(edge.network))
  g.dist <- (g.dist[enm, 1:max.node.id] - 1) / 2
  if(lookup.mat) {
    # 3d lookup matrix (edge node 1, edge node 2, focal node)
    m <- array(NA, rep(vcount(g), 3))
    g.el <- get.edgelist(g, names = F)
    # this creates an index for the 3d matrix and merges it with the melted
    b <- do.call(rbind, replicate(ncol(g.dist), g.el, simplify=F)) %>%
      cbind(reshape2::melt(g.dist)) %>%
      select(-Var1)
    m[b %>% select(`1`, `2`, `Var2`) %>% as.matrix] <- b$value
    return(m)
  }
  if(add.names) {
    if("matrix" %in% class(g.dist)) {
      if(is.null(ns)) {
        ns <- 1:ncol(g.dist)
      }
      colnames(g.dist) <- ns
    } else {
      names(g.dist) <- ns
    }
    g.el <- get.edgelist(g)
    rownames(g.dist) <- paste(g.el[,1],g.el[,2],sep="-")
  }
  return(g.dist)
}



