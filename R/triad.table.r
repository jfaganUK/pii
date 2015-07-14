#' Triad Table
#'
#' Creates a table that contains triadID - the three nodes that make up the triad, nodeID - the character name of the node,
#' direction of the triad from the node, the valence of the triad from the node's perspective, distance of the triad from the node,
#' and nodeNum - the numderic value of the node
#'
#' @param g An igraph graph
#' @export


triad.table <- function(g) {
  triad_table <-
    data.table(
      triadID = character(), nodeID = character(), direction = character(), valence = integer(),
      distance = numeric(), nodeNum = integer()
    )
  triad_holder = cliques(g, min = 3, max = 3)
  e.dist <- edge.distance(g)
  node_distance_holder <- shortest.paths(g, V(g))
  for (i in 1:length(triad_holder)) {
    cat(' Triad: ', i, '\n')
    node1 = V(g)[triad_holder[[i]][1]]$name
    node2 = V(g)[triad_holder[[i]][2]]$name
    node3 = V(g)[triad_holder[[i]][3]]$name
    triname = paste(node1, node2, node3, sep = "-")
    edge1_2 = paste(node1, node2, sep = "-")
    edge1_3 = paste(node1, node3, sep = "-")
    edge2_3 = paste(node2, node3, sep = "-")
    for (number in 1:vcount(g)) {
      refnode = V(g)[number]$name
      dis1 = e.dist[edge1_2, refnode]
      dis2 = e.dist[edge1_3, refnode]
      dis3 = e.dist[edge2_3, refnode]

      # inward facing triad
      if (dis1 == dis2 & dis1 == dis3) {
        nodedis1_2 = node_distance_holder[node1, node2]
        nodedis1_3 = node_distance_holder[node1, node3]
        nodedis2_3 = node_distance_holder[node2, node3]
        if (nodedis1_2 !=  nodedis1_3) {
          if (nodedis1_2 != nodedis2_3) {
            ref = edge1_2
          }else{
            ref = edge1_3
          }
        }else{
          ref = edge2_3
        }
        if (ref == edge1_2) {
          truev = valence_set[1]
        }
        if (ref == edge1_3) {
          truev = valence_set[3]
        }
        if (ref == edge2_3) {
          truev = valence_set[2]
        }
        row <-
          data.table(
            triadID = triname, nodeID = refnode, direction = 'IN ', valence = truev,
            distance = e.dist[edge1_2, refnode], nodeNum = as.integer(V(g)[refnode])
          )
        triad_table <- rbind(triad_table, row, use.names = T)
      }

      # outward facing triad
      else{
        if (dis1 != dis2) {
          if (dis1 != dis3) {
            diff = edge1_2
          }else{
            diff = edge1_3
          }
        }else{
          diff = edge2_3
        }
        p <- c(triad_holder[[i]], triad_holder[[i]][1])
        valence_set <- E(g, path = p)$valence
        #print(valence_set)
        if (diff == edge1_2) {
          truev = valence_set[1]
        }
        if (diff == edge1_3) {
          truev = valence_set[3]
        }
        if (diff == edge2_3) {
          truev = valence_set[2]
        }
        row <-
          data.table(
            triadID = triname, nodeID = refnode, direction = 'OUT', valence = truev,
            distance = e.dist[diff, refnode],  nodeNum = as.integer(V(g)[refnode])
          )
        triad_table <- rbind(triad_table, row, use.names = T)
      }
    }
  }
  return(triad_table)
}
