crosses<- function(g, pii.beta = -0.7, tolerance = 0.000001){
  n1 <- c()
  n2 <- c()
  table <- c()

  # get the distance vectors
  e.dist <- edge.distance(g)
  e.dist <- matrix(as.integer(e.dist), nrow=nrow(e.dist)) # convert to an integer matrix
  max.distance <- max(e.dist)
  edgevalence <- E(g)$valence
  vec <- valenceEdgeCounts(e.dist, edgevalence, max.distance)

  max.degree <- max(degree(g, mode='total'))
  for(i in 1:(vcount(g)-1)){
    for(j in (i+1):vcount(g)){
      if(i != j){
        beta.vec1 <- overlaps(vec, pii.beta, max.degree, i, j)

        if(length(beta.vec1) > 0) {
          cat('i: ', i, 'j: ', j, '\n')
        }
        o <- numeric(length(beta.vec1))
        ## Because pii.x has beta in it, we use an initial value of beta to start and get an estimate
        ## of the beta overlap. Then we use the values of the estimated beta overlaps to calculate
        ## the positions again.
        if(length(beta.vec1) > 0) {
           # output vector
          for(k in 1:length(beta.vec1)) {
            prev.beta <- pii.beta
            curr.beta <- mean(c(prev.beta, beta.vec1[k])) # find a mid-point between
            while(abs(curr.beta - prev.beta) > tolerance) {
              prev.beta <- curr.beta
              (beta.vec <- overlaps(vec, curr.beta, max.degree, i, j))
              if(length(beta.vec) > 1) {
                curr.beta <- mean(c(prev.beta, beta.vec[k]))
              }
              if(is.na(curr.beta)){
                break
              }
            }
            o[k] <- curr.beta
          }
        }



#         xx <- polyroot((vec$pos[,i] - vec$neg[,i]) - (vec$pos[,j] - vec$neg[,j]))
         length.table <- length(table)
#         table <- c(table, Re(xx[round(Im(xx), digits = 12) == 0]))
#         table <- table[table < -0.5 & table > -0.9]
         table <- c(table, o)
        #while(curr.beta-prev.beta > tolerance){
          #prev.beta <- curr.beta
          #curr.beta <- polyroot((vec$pos[,i] - vec$neg[,i]) - (vec$pos[,j] - vec$neg[,j]))
        #}
         length.curr.table <- length(table)
         if(length.table != length.curr.table){
           n1 <- c(n1, rep(i, (length.curr.table - length.table)))
           n2 <- c(n2, rep(j, (length.curr.table - length.table)))
         }
      }
    }
  }
  na.omoit(data_frame("n1" = n1, "n2" = n2, "beta" = table,
             "graphid" = ifelse(is.null(get.graph.attribute(g, "graphid")), "", get.graph.attribute(g, "graphid"))))
}


overlaps <- function(vec, beta, max.degree, i, j){
  pii.x <- (log(2) - log(abs(beta))) / log(max.degree)
  xx <- polyroot((vec$pos[,i] ^ pii.x - vec$neg[,i] ^ pii.x) - (vec$pos[,j]^pii.x - vec$neg[,j]^pii.x))
  Re(xx[which(Re(xx) < 0)[1]])
  #tt <- Re(xx[round(Im(xx), digits = 12) == 0])
  #tt[tt <= -0.5 & tt >= -0.9]
  #Re(xx)[1]
}
