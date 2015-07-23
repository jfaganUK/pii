library(igraph)
library(reshape2)
library(data.table)
library(ggplot2)
library(Rcpp)
library(rbenchmark)
library(scales)
library(pii)
library(dplyr)
library(magrittr)
library(gridExtra)
library(parallel)

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

#g <- graph.edgelist(as.matrix(samp.el[,c('source','target'),with=F]), directed = T)
g <- graph.data.frame(samp.el)
E(g)$valence <- samp.el$valence
E(g)$weight <- samp.el$value
g <- as.undirected(g, edge.attr.comb = 'min')
#g <- simplify(g, edge.attr.comb = 'min')
g.samp <- g

countTieValence <- function(v, val) {
  sum(E(g.samp)[inc(v)]$valence == val)
}

g.samp <- remove.edge.attribute(g.samp, "value")
save(g.samp, file = 'g.samp.rda')

tt <- data.table(nm = V(g.samp)$name,
           pii = pii(g.samp),
           pos.ties = sapply(V(g.samp), function(v) { countTieValence(v, 1) }),
           neg.ties = sapply(V(g.samp), function(v) { countTieValence(v, -1) }))
tt[, ratio := pos.ties / (pos.ties + neg.ties)]
tt[order(ratio, decreasing=T)]

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
  el <- matrix(c('a','b','a','c','b','c', 'a','d', 'd','e', 'd','f', 'd','g', 'g','e'), ncol=2, byrow=T)
  g <- graph.edgelist(el, directed = F)
  E(g)$valence <- 1
  E(g)[1 %--% 2]$valence <- -1
  E(g)[1 %--% 3]$valence <- -1
  E(g)[5 %--% 7]$valence <- -1
  E(g)$color <- ifelse(E(g)$valence == -1, "red", "black")
  g
}

getThinkGraph2 <- function() {
  el <- matrix(c('a','b', 'a','c', 'b','c', 'b','d', 'c','d', 'd','e', 'e','f'), ncol=2, byrow=T)
  g <- graph.edgelist(el, directed = F)
  E(g)$valence <- 1
  E(g)[5 %--% 6]$valence <- -1
  E(g)$color <- ifelse(E(g)$valence == -1, "red", "black")
  g
}
g <- getThinkGraph2()
plot(g)
pii(g, triadic=T)

cliques(g, min = 3, max = 3)
triad_holder = cliques(g, min = 3, max = 3)
node_distance_holder <- shortest.paths(g, V(g))





# piis <- data.table(node = character(), pii = numeric(), delta = numeric())
# for(d in seq(0,1,by=0.1)) {
#   cat('Delta: ', d, '\n')
#   td <- triadicPII(g, pii.delta = d)
#   piis <- rbind(piis, data.table(node = V(g)$name, pii = td, delta = d))
# }

ggplot(piis) + geom_line(aes(x=delta, y=pii, group=node, color=node))



g1 <- randomGraph()
V(g1)$name <- paste0(V(g1)$name, '1')
V(g1)$name <- paste0(V(g1)$name, '1')
g2 <- randomGraph()
V(g2)$name <- paste0(V(g2)$name, '2')
igraph.options(edge.attr.comb = list(valence = max))
g3 <- graph.union(g1, g2)

E(g3)$color <- ifelse(!is.na(E(g3)$color_1), E(g3)$color_1, E(g3)$color_2)
plot(g3)





piis <- data.table(nd=character(), pii=numeric(), beta=numeric())
for(b in seq(-1, 0, by=0.001)) {
  piis <- rbind(piis, data.table(nd = V(g)$name, pii=pii(g, pii.beta=b), beta = b))
}
ggplot(piis) + geom_line()


#centralization.degree(g)$centralization
#vcount(g)
#ecount(g)
#graph.density(g)
#length(E(g)[valence == -1])
#length(E(g)[valence == 1])
#average.path.length(g)
###var(pii)
#beta.sequence <- seq(-1, -0.1, by=0.1)
#all.pii <- data.table(node=character(), pii.value = numeric(), beta = numeric())
#g <- all.ring.graphs[[1]]
#g.ed <- edge.distance(g)
#for(b in beta.sequence) {
#  g.pii <- pii(g,e.dist = g.ed, pii.beta = b)
#  td <- data.table(node=1:vcount(g), pii.value=as.numeric(g.pii), degree = degree(g), beta=b)
#  all.pii <- rbind(all.pii, td)
#}


#md(all.ring.graphs[[1]])
gp <- do.call('rbind', mclapply(all.graphs, graphData, mc.cores=5))
ggplot(gp, aes(y=rankCor, x=diameter)) + geom_point() + geom_smooth(method='lm')

### Random graphs - Watts Strogatz ############################################

g <- randomGraph()
plot(g)

rand.graphs = list()
for(i in 1:10){
  rand.graphs[[i]] <- randomGraph(i)
}

gp <- do.call('rbind', mclapply(rand.graphs, graphData, mc.cores=5))
ggplot(gp, aes(y=rankCor, x=lowCorBeta)) + geom_point() + geom_smooth(method='lm')
ggplot(gp, aes(y=rankCor, x=propNegEdge)) + geom_point() + geom_smooth(method='lm')
summary(lm1 <- lm(rankCor ~ meanTrans + avgPathLength + degCentralization + propNegEdge + density + modularity, data=gp))
summary(lm2 <- lm(rankCor ~ avgPathLength + propNegEdge + avgPathLength:numNegEdge, data=gp))
summary(lm2 <- lm(rankCor ~ modularity, data=gp))




library(parallel)
all.graphs <- list()
for(i in 3:10) {
  all.graphs <- c(all.graphs, getAllRingGraphs(i))
  all.graphs <- c(all.graphs, getAllRingGraphs(i, chain=T))
}
nd <- do.call('rbind', mclapply(all.graphs, nodeData, mc.cores=5))

# cnd <- do.call('rbind', lapply(getAllRingGraphs(3, chain = T), nodeData))
# for(i in 4:10){
#   cnd <- rbind(cnd, do.call('rbind', lapply(getAllRingGraphs(i, chain = T), nodeData)))
# }


#' table 3
#' graph level
#'


nd3 <- group_by(nd, graphid, beta) %>%
  summarize(meanPII = mean(pii.value), varPII = var(pii.value))


all.graph.name <- unique(nd$graphid)
k <- 31
pdf(paste0('catalog.pdf'))
for(k in 1:length(all.graph.name))  {
  plot(all.graphs[[all.graph.name[k]]])
  p <- ggplot(nd[graphid == all.graph.name[k]]) +
    geom_line(size=3, alpha=0.6, position=position_jitter(height=0.02, width=0.02),
              aes(x=beta, y=pii.value, group=nodeid, color=nodeid))
  print(p)
}
dev.off()


#centralization.degree(g)$centralization
#vcount(g)
#ecount(g)
#graph.density(g)
#length(E(g)[valence == -1])
#length(E(g)[valence == 1])
#average.path.length(g)
###var(pii)
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

all.pii[, is.neg := ifelse(pii.value < 0, "neg", "pos")]
all.pii[, pii.value := rescale(pii.value)]
ggplot(all.pii, aes(x=beta, y=delta, fill=pii.value)) + geom_tile() + facet_wrap(~ node) +
  scale_x_continuous('Beta') + scale_y_continuous('Delta') +
  scale_fill_gradient2()

all.pii.agg <- all.pii[,list(m.pii = mean(pii.value, na.rm=T)), by=list(beta, delta)]
ggplot(all.pii.agg, aes(x=beta, y=delta, fill=m.pii)) + geom_tile() +
  scale_x_continuous('Beta') + scale_y_continuous('Delta') +
  scale_fill_gradient2(low = "white", mid="brown", high = "steelblue")

ggplot(all.pii.agg) + geom_point(size=15, aes(x = beta, y=delta, color=m.pii))


library(parallel)
g <- randomGraph()
plot(g)

tb <- -1
p1 <- pii(g, pii.beta = tb)
x <- do.call('rbind', mclapply(seq(-1, -0.01, by=0.01), function(b) {
  p <- pii(g, pii.beta = b)
  data.table(b = b, rc = cor(p1, p, method = "spearman"))
}, mc.cores = 6))
ggplot(x, aes(x=b, y=rc)) + geom_line() + geom_smooth()

#run calc, find beta where rc goes below .707
#calculate avg minimum distance to negative edge
#calculate avg distance to negative edge

#finds beta where rc goes below .707
