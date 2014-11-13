pii
===

R Package for the Political Independence Index.


```R
# Intall PII in R
installed.packages("devtools")
devtools::install_github("jfaganUK/pii")

# Load the package
library(pii)
# The pii package requires igraph. The network object must be an igraph object
library(igraph)

# create some random graph
N <- 10
g <- erdos.renyi.game(N, runif(N))

# add some random valences to each edge
# valence can be 1 (positive tie) or -1 (negative tie)
# the 'valence' edge attribute is required
E(g)$valence <- sample(c(-1, 1), length(E(g)), replace=T)
V(g)$name <- letters[1:N]

pii(g)
```
