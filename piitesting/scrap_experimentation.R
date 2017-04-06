triad_table
triadsPerNode <- c(1:length(unique(triad_table$focalNode)))
numTriadPerNode <- lapply(triadsPerNode, function(n){nrow(triad_table[triad_table$focalNode == n, ])})
numTriadPerNode <- as.matrix(numTriadPerNode)
