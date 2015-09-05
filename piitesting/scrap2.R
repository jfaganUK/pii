e.dist <- edge.distance(g)

e.dist.l <- list(ncol(e.dist))
for(j in 1:ncol(e.dist)) {
  m1 <- matrix(-1, nrow = vcount(g), ncol = vcount(g))
  rownames(m1) <- colnames(e.dist)
  colnames(m1) <- colnames(e.dist)

  nms <- strsplit(rownames(e.dist), "-")
  for(i in 1:nrow(e.dist)) {
    x <- nms[[i]]
    m1[x[1], x[2]] <- e.dist[i, j]
    m1[x[2], x[1]] <- e.dist[i, j]
  }

  e.dist.l[[j]] <- m1
}

cppFunction('int getEdgeDistance(List edist, int n1, int n2, int j) {
NumericMatrix m = edist[j];
return m(n1, n2);
            }')

getEdgeDistance(e.dist.l, 2, 5, 26)
