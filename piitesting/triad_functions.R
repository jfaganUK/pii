getTriadRow <- function(i) {
  refnode = V(g)[i]$name
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
        distance = e.dist[edge1_2, refnode], nodeNum = as.integer(V(g)[refnode]))
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
    truev <- edge.table[diff, valence]
    row <-
      data.table(
        triadID = triname, nodeID = refnode, direction = 'OUT', valence = truev,
        distance = e.dist[diff, refnode],  nodeNum = as.integer(V(g)[refnode])
      )
  }
  return(row)
}

do.call('rbind', lapply(1:vcount(g), 'getTriadRow'))



getTriadRows <- function(j) {
  cat(' Triad: ', j, '\n')
  node1 = V(g)[triad_holder[[j]][1]]$name
  node2 = V(g)[triad_holder[[j]][2]]$name
  node3 = V(g)[triad_holder[[j]][3]]$name
  triname = paste(node1, node2, node3, sep = "-")
  edge1_2 = paste(node1, node2, sep = "-")
  edge1_3 = paste(node1, node3, sep = "-")
  edge2_3 = paste(node2, node3, sep = "-")

  rows <- do.call('rbind', lapply(1:vcount(g), 'getTriadRow'))
  return(rows)
}
