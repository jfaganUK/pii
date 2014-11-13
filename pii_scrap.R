library(igraph)
library(reshape2)
library(data.table)
library(ggplot2)
library(Rcpp)
library(rbenchmark)
library(scales)

options(stringsAsFactors=F)

###### trim ##############################
# trims the leading and trailing spaces from a string (you can provide a vector)
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

plotty <- function(g) {
  if(!('igraph' %in% class(g))) {
    stop('Requires an igraph object.')
  } else {
    igraph::plot.igraph(g, vertex.size=8, vertex.color='SkyBlue2', vertex.frame.color='SkyBlue2', edge.arrow.size=0.5,
         vertex.label.family='Open Sans', vertex.label.cex=0.5)
  }
}

# create some random graph
N <- 10
g <- erdos.renyi.game(N, runif(N))

# add some random valences to each edge
# valence can be 1 (positive tie) or -1 (negative tie)
# the 'valence' edge attribute is required
#E(g)$valence <- sample(c(-1, 1), length(E(g)), replace=T)
V(g)$name <- letters[1:N]

pii(g)


# there are 18 nodes
x <- read.csv('sampson.txt',stringsAsFactors=F)
N <- 18
networkNames <- c("samplk1","samplk2","samplk3","sampdlk","sampes","sampdes","sampin","sampnin","samppr","sampnpr")
x$Node <- trim(x$Node)

samp <- list()
for(i in 1:length(networkNames)) {
  y <- x[(1:18) + ((i - 1) * 18),]
  rownames(y) <- y$Node
  y <- as.matrix(y[,-1])
  samp[[i]] <- y
}
names(samp) <- networkNames


g <- igraph::graph.adjacency(samp[['samplk3']])
plotty(g)

# indegree for everything...
do.call('rbind', lapply(samp, function(x) {
  g <- graph.adjacency(x)
  degree(g, mode="in")
}))

# The political networks paper uses Esteem and Disesteem (sampes and sampdes)

g.pos <- graph.adjacency(samp[['sampes']])
g.neg <- graph.adjacency(samp[['sampdes']])

samp.el <- rbind(data.table(melt(samp[['sampes']]))[value != 0,][,valence := 1],
                 data.table(melt(samp[['sampdes']]))[value != 0,][,valence :=-1])
setnames(samp.el, c('Var1','Var2'),c('source','target'))



g <- graph.edgelist(as.matrix(samp.el[,c('source','target'),with=F]), directed = F)
E(g)$valence <- samp.el$valence
E(g)$weight <- samp.el$value
g <- simplify(g, edge.attr.comb = 'min')
g.samp <- g

g.pos <- graph.edgelist(as.matrix(samp.el[valence > 0,c('source','target'),with=F]), directed = F)
g.neg <- graph.edgelist(as.matrix(samp.el[valence < 0,c('source','target'),with=F]), directed = F)
E(graph.intersection(g.pos, g.neg))


### Testing PII ###############################################################

### The small example from the paper
# g <- as.undirected(graph.empty()) + vertices(letters[1:9])
# g <- g + path('i','e','a','c','g')
# g <- g + path('h','d','a','b','f')
# E(g)$valence <- 1
# E(g)[vid('a') %--% vid('c')]$valence <- -1
# E(g)[vid('a') %--% vid('b')]$valence <- -1

### An example (figure 2) from the paper
getThinkGraph <- function() {
  el <- matrix(c('a','b','a','c','b','c', 'a','d', 'd','e', 'd','f', 'd','g'), ncol=2, byrow=T)
  g <- graph.edgelist(el, directed = F)
  E(g)$valence <- 1
  E(g)[1 %--% 2]$valence <- -1
  E(g)[1 %--% 3]$valence <- -1
  g
}
g <- getThinkGraph()

t1 <- proc.time()
beta.sequence <- seq(-1, -0.1, by=0.1)
delta.sequence <- seq(0, 1, by=0.1)
all.pii <- data.table(node=character(), pii.value = numeric(), beta = numeric(), delta = numeric())
g.ed <- edge.distance(g)
for(i in 1:length(beta.sequence)) {
  for(j in 1:length(delta.sequence)) {
    cat('beta: ', i, ' - delta: ', j,'\r')
    samp.pii <- pii(g,e.dist = g.ed, pii.beta = beta.sequence[i], triadic = T, pii.delta = delta.sequence[j])
    all.pii <- rbind(all.pii,
                     data.table(node=names(samp.pii), pii.value=as.numeric(samp.pii),
                                beta=beta.sequence[i], delta=delta.sequence[j]))
  }
}
proc.time() - t1
all.pii[, is.neg := ifelse(pii.value < 0, "neg", "pos")]
all.pii[, pii.value := rescale(pii.value)]
ggplot(all.pii, aes(x=beta, y=delta, fill=pii.value)) + geom_tile() + facet_wrap(~ node) +
  scale_x_continuous('Beta') + scale_y_continuous('Delta') +
  scale_fill_gradient(low = "white", high = "steelblue")

all.pii.agg <- all.pii[,list(m.pii = mean(pii.value, na.rm=T)), by=list(beta, delta)]
ggplot(all.pii.agg, aes(x=beta, y=delta, fill=m.pii)) + geom_tile() +
  scale_x_continuous('Beta') + scale_y_continuous('Delta') +
  scale_fill_gradient(low = "white", high = "steelblue")

ggplot(all.pii.agg) + geom_point(size=15, aes(x = beta, y=delta, color=m.pii))


### Compute an auto-correlation matrix

d <- 0.5
all.pii.d <- all.pii[delta == d]
psi.corr <- data.table(beta1 = numeric(), beta2 = numeric(), psi = numeric())
for(b1 in unique(all.pii.d$beta)) {
  for(b2 in unique(all.pii.d$beta)) {
    psi <- cor(all.pii.d[beta == b1]$pii.value, all.pii.d[beta == b2]$pii.value)
    psi.corr <- rbind(psi.corr, data.table(beta1=b1, beta2=b2, psi=psi))
  }
}

ggplot(psi.corr, aes(x=beta1, y=beta2, fill=psi)) + geom_tile()

psi.mat <- matrix(psi.corr[order(beta1, beta2)]$psi, nrow=10)

eigen(psi.mat)


### benchmark different versions
e.dist <- edge.distance(g)
benchmark(pii(g), pii2(g), pii(g, e.dist = e.dist), pii2(g, e.dist = e.dist), pii2(g, e.dist = e.dist, triadic = T), replications=10, order='relative')
benchmark(edge.distance(g), edge.distance2(g), replications = 10, order = 'relative')


g.lo <- layout.kamada.kawai(g)

x <- pii2(g, pii.beta = -0.8)
plot(g, layout=g.lo, vertex.size=rescale(x)*50, vertex.color='SkyBlue2',
     vertex.frame.color='SkyBlue2', edge.arrow.size=0.5,
     vertex.label.family='Open Sans', vertex.label.cex=2)



