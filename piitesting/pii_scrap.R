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
#Re(xx[round(Im(xx), digits = 12) == 0])
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

data('g.samp')

countTieValence <- function(v, val) {
  sum(E(g.samp)[inc(v)]$valence == val)
}

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
# ggplot(piis) + geom_line(aes(x=delta, y=pii, group=node, color=node))

# g1 <- randomGraph()
# V(g1)$name <- paste0(V(g1)$name, '1')
# V(g1)$name <- paste0(V(g1)$name, '1')
# g2 <- randomGraph()
# V(g2)$name <- paste0(V(g2)$name, '2')
# igraph.options(edge.attr.comb = list(valence = max))
# g3 <- graph.union(g1, g2)
#
# E(g3)$color <- ifelse(!is.na(E(g3)$color_1), E(g3)$color_1, E(g3)$color_2)
# plot(g3)


piis <- data.table(nd=character(), pii=numeric(), beta=numeric())
for(b in seq(-1, 0, by=0.001)) {
  piis <- rbind(piis, data.table(nd = V(g)$name, pii=pii(g, pii.beta=b), beta = b))
}
ggplot(piis, aes(x=beta, y=pii, group=nd, color=nd)) + geom_line()




#centralization.degree(g)$centralization
#vcount(g)
#ecount(g)
#graph.density(g)
#length(E(g)[valence == -1])
#length(E(g)[valence == 1])
#average.path.length(g)
###var(pii)
# beta.sequence <- seq(-1, -0.1, by=0.1)
# all.pii <- data.table(node=character(), pii.value = numeric(), beta = numeric())
# g <- all.ring.graphs[[1]]
# g.ed <- edge.distance(g)
# for(b in beta.sequence) {
#  g.pii <- pii(g,e.dist = g.ed, pii.beta = b)
#  td <- data.table(node=1:vcount(g), pii.value=as.numeric(g.pii), degree = degree(g), beta=b)
#  all.pii <- rbind(all.pii, td)
# }
g <- getThinkGraph2()
beta.sequence <- seq(-1, -0.1, by=0.01)
all.pii <- data.table(node=character(), pii.value = numeric(), beta = numeric())
g.ed <- edge.distance(g)
for(b in beta.sequence) {
  g.pii <- pii(g,e.dist = g.ed, pii.beta = b)
  td <- data.table(node=V(g)$name, pii.value=as.numeric(g.pii), beta=b)
  all.pii <- rbind(all.pii, td)
}

ggplot(all.pii,aes(x=beta, y=pii.value, group = node, color = node)) +
  scale_color_manual(values = c('a' = '#E41A1C', 'b' = '#377EB8', 'c' = '#4DAF4A', 'd' = '#984EA3', 'e' = '#FF7F00', 'f' = '#FFFF33')) +
  geom_line(alpha = 0.8, size = 2) + thm

V(g)$color <- RColorBrewer::brewer.pal(6, 'Set1')

### Animation the changes
library(animation)
saveGIF({
  par(mar = c(0,0,2.5,0), bg='#eeeeee' )
  lo <- layout_nicely(g)
  for(b in c(beta.sequence, rev(beta.sequence))) {
    V(g)$size <- all.pii %>% filter(beta == b) %>%
      select(pii.value) %>% collect %>% .[[1]] %>% rescale(to = c(5,35))
    plot(g, layout = lo)
    title(paste0('Beta = ', b))
  }
  dev.off()
}, interval = 0.01, movie.name = 'pii_by_beta.gif', ani.width = 600, ani.height = 600)

#md(all.ring.graphs[[1]])
# gp <- do.call('rbind', mclapply(all.graphs, graphData, mc.cores=5))
# ggplot(gp, aes(y=rankCor, x=diameter)) + geom_point() + geom_smooth(method='lm')

### Random graphs - Watts Strogatz ############################################

par(mar=c(0,0,0,0))
g <- randomGraph()
plot(g)

rand.graphs = list()
for(i in 1:10){
  rand.graphs[[i]] <- randomGraph()
}

gp <- do.call('rbind', mclapply(rand.graphs, graphData, mc.cores=5))
save(gp, file = "./piitesting/gp.rda")

load('./piitesting/gp.rda')
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
ggplot(x, aes(x=b, y=rc)) + geom_line() + geom_smooth(method='loess') +
  scale_y_continuous(lim=c(0,1.0))


### Two networks ###############################################################
# one that's stable, and one that isn't

# An unstable graph
g1 <- getThinkGraph2()
par(mar=c(0,0,0,0))
plot(g1)
pii.beta1 <- optimal.rank.beta(g1)
suppressWarnings(cor(pii(g1, pii.beta = -0.9), pii(g1, pii.beta=-0.5), method = "spearman"))
pii.diagnostic.plot(g1)

# A stable graph
g2 <- g1
E(g2)$valence <- 1
E(g2)[1 %--% 2]$valence <- -1
E(g2)$color <- ifelse(E(g2)$valence == -1, "red", "black")
plot(g2)
(pii.beta2 <- optimal.rank.beta(g2, comp.left = -0.9, comp.right = -0.5))
suppressWarnings(cor(pii(g2, pii.beta = -0.9), pii(g2, pii.beta=-0.5), method = "spearman"))
pii.diagnostic.plot(g2)




### Beta - Node Stability Graph ################################################
g <- randomGraph(maxnegtie = 0.2, badeggchance = 0.1)
gp <- graphData(g)
par(bg = '#eeeeee', mar = c(0,0,0,0))
plot(g, vertex.size = 4, vertex.label = NA, edge.width = 1, vertex.color = '#333333')
dev.off(); dev.new()

inc <- 0.01
x <- do.call('rbind', mclapply(seq(-1, -0.01, by=inc), function(b) {
  p <- pii(g, pii.beta = b)
  piir <- rank(p)
  data.table(b = b, nd = V(g)$name, pii = p, piir = piir)
}, mc.cores = 6))

ggplot(x, aes(x=b, y=pii, group=nd, color=nd)) +
  scale_x_continuous(expression(beta), limits = c(-1,-0.5)) +
  scale_y_continuous('PII Score', limits = c(-10, 5)) +
  geom_line(size=1) + thm + theme(legend.position = 'none')

ggplot(x, aes(x=b, y=piir, group=nd, color=nd)) +
  geom_line(alpha=0.9, size = 1) +
  scale_x_continuous(expression(beta), lim = c(-1, -0.5)) +
  scale_y_continuous('PII Ranked Score') +
  thm + theme(legend.position = 'none')

comp.left <- -0.9
comp.right <- -0.5
rc <- do.call('rbind', lapply(unique(x$b)[-1], function(bb) {
  pii <- x[b == bb, pii]
  pii.left <- x[b == comp.left, pii]
  pii.right <- x[b == comp.right, pii]
  data.table(b = bb, rc.left = cor(pii, pii.left, method = "spearman"), rc.right = cor(pii, pii.right, method = "spearman"))
}))

rc.diff <- function(b, g, comp.left=-0.9, comp.right=-0.5) {
  p <- pii(g, pii.beta = b)
  pii.left <- pii(g, pii.beta = comp.left)
  pii.right <- pii(g, pii.beta = comp.right)
  rc.left = cor(p, pii.left, method = "spearman")
  rc.right = cor(p, pii.right, method = "spearman")
  return(abs(rc.left - rc.right))
}
rc.diff(-0.63, g)
o <- optim(par=c(-0.8), fn = rc.diff, gr = NULL,
           g = g, comp.right = comp.right, comp.left = comp.left,
           method='Brent', lower = -1, upper = -0.001)
cross.y <- cor(pii(g, pii.beta = o$par), pii(g, pii.beta = comp.left), method = "spearman")
cross.point <- data.frame(xx=o$par, yy=cross.y)

## Need to find the crossing
rc.m <- melt(rc, id.vars='b')
ggplot() +
  geom_line(data = rc.m, aes(x=b, y=value, group=variable)) +
  geom_text(data=cross.point, aes(x=xx, y=yy, label = round(xx,3)), size = 4, vjust=-2) +
  scale_y_continuous('Rank Correlation', lim=c(0,1.0)) +
  scale_x_continuous(expression(beta), lim=c(-1, -0.01)) + thm

  ### Node Crossing ##############################################################
# g <- randomGraph()
g <- getThinkGraph()
plot(g)

nm <- V(g)$name %>% combn(2) %>% t %>% data.table %>%
  rename(name1 = V1, name2 = V2)

i <- 11
# nm1 <- nm$name1[i]
# nm2 <- nm$name2[i]
nm1 <- 'f'
nm2 <- 'g'

pii.diff <- function(b, nm1, nm2) {
  p <- pii(g, pii.beta = b)
  # p <- rank(p)
  (p[nm1] - p[nm2])^2
}
po <- optim(par = -0.1, pii.diff, gr = NULL, nm1 = nm1, nm2 = nm2,
            method = 'Brent', lower = -0.999, upper = -0.001)
p <- pii(g, pii.beta = po$par)
# p <- rank(p)
p[nm1]
p[nm2]


### All graphs analysis ########################################################

all.graphs <- mclapply(1:300, function(i) {
  return(randomGraph(i, maxnegtie = 0.1))
}, mc.cores=5)
gd <- do.call('rbind', mclapply(all.graphs, graphData, mc.cores=5))


ggplot(gd) +
  geom_bar(aes(x=optimBeta), binwidth=0.005)

ggplot(gd, aes(x=rankCor, y=optimBeta)) +
  geom_point() + geom_smooth(method='lm')

lm(optimBeta ~ avgMinDistToNegEdge + avgPathLength + modularity, data = gd) %>%
  summary()

select(gd, -graphid) %>% cor %>% melt %>%
  filter(Var1 == 'optimBeta' & value < 1 ) %>% arrange(value) %>%
  mutate(value = round(value, 2))

library(randomForest)
rf <- randomForest(optimBeta ~ diameter + density + meanDegree + nodeCount +
                     edgeCount + avgPathLength + propNegEdge + modularity +
                     meanTrans + avgMinDistToNegEdge + avdDistOfNegEdge,
                   data = na.omit(gd))

varImpPlot(rf)


#finds beta where rc goes below .707

all.graphs <- list()
for(i in 3:10) {
  all.graphs <- c(all.graphs, getAllRingGraphs(i))
  all.graphs <- c(all.graphs, getAllRingGraphs(i, chain=T))
}
nd <- do.call('rbind', mclapply(all.graphs, nodeData, mc.cores=5))

ring.graphs = list()
for(i in 3:10){
  ring.graphs <- c(ring.graphs, getAllRingGraphs(num = i))
}
rd <- do.call('rbind', mclapply(ring.graphs, graphData, mc.cores=5))
save(rd, file = "./piitesting/all-ring-03-10-stats.rda")


chain.graphs = list()
for(i in 3:10){
  chain.graphs <- c(chain.graphs, getAllRingGraphs(num = i, chain = T))
}
cd <- do.call('rbind', mclapply(chain.graphs, graphData, mc.cores=5))
save(cd, file = "./piitesting/all-chain-03-10-stats.rda")


rand.graphs = list()
for(i in 1:300){
  rand.graphs[[i]] <- randomGraph(maxnegtie = 0, pendchance = 0.00, badeggchance = 0.00)
}

rga <- do.call('rbind', mclapply(rand.graphs, graphData, mc.cores=5))
save(rga, file = "./piitesting/random-graph-a-stats.rda")


rand.graphs = list()
for(i in 1:300){
  rand.graphs[[i]] <- randomGraph(maxnegtie = 0.05, pendchance = 0.00, badeggchance = 0.00)
}

rg_b <- do.call('rbind', mclapply(rand.graphs, graphData, mc.cores=5))
save(rg_b, file = "./piitesting/random-graph-b-stats.rda")


rand.graphs = list()
for(i in 1:300){
  rand.graphs[[i]] <- randomGraph(maxnegtie = 0.1, pendchance = 0.00, badeggchance = 0.00)
}

rgc <- do.call('rbind', mclapply(rand.graphs, graphData, mc.cores=5))
save(rgc, file = "./piitesting/random-graph-c-stats.rda")

rand.graphs = list()
for(i in 1:300){
  rand.graphs[[i]] <- randomGraph(maxnegtie = 0.2, pendchance = 0.00, badeggchance = 0.00)
}

rgd <- do.call('rbind', mclapply(rand.graphs, graphData, mc.cores=5))
save(rgd, file = "./piitesting/random-graph-d-stats.rda")

rand.graphs = list()
for(i in 1:300){
  rand.graphs[[i]] <- randomGraph(maxnegtie = 0.1, pendchance = 0.01, badeggchance = 0.00)
}

rge <- do.call('rbind', mclapply(rand.graphs, graphData, mc.cores=5))
save(rge, file = "./piitesting/random-graph-e-stats.rda")

rand.graphs = list()
for(i in 1:300){
  rand.graphs[[i]] <- randomGraph(maxnegtie = 0.1, pendchance = 0.02, badeggchance = 0.00)
}

rgf <- do.call('rbind', mclapply(rand.graphs, graphData, mc.cores=5))
save(rgf, file = "./piitesting/random-graph-f-stats.rda")

rand.graphs = list()
for(i in 1:300){
  rand.graphs[[i]] <- randomGraph(maxnegtie = 0.1, pendchance = 0.05, badeggchance = 0.00)
}

rgg <- do.call('rbind', mclapply(rand.graphs, graphData, mc.cores=5))
save(rgg, file = "./piitesting/random-graph-g-stats.rda")

rand.graphs = list()
for(i in 1:300){
  rand.graphs[[i]] <- randomGraph(maxnegtie = 0.1, pendchance = 0.1, badeggchance = 0.00)
}

rgh <- do.call('rbind', mclapply(rand.graphs, graphData, mc.cores=5))
save(rgh, file = "./piitesting/random-graph-h-stats.rda")

rand.graphs = list()
for(i in 1:300){
  rand.graphs[[i]] <- randomGraph(maxnegtie = 0.1, pendchance = 0.2, badeggchance = 0.00)
}

rgi <- do.call('rbind', mclapply(rand.graphs, graphData, mc.cores=5))
save(rgi, file = "./piitesting/random-graph-i-stats.rda")

rand.graphs = list()
for(i in 1:300){
  rand.graphs[[i]] <- randomGraph(maxnegtie = 0.1, pendchance = 0.02, badeggchance = 0.01)
}

rgj <- do.call('rbind', mclapply(rand.graphs, graphData, mc.cores=5))
save(rgj, file = "./piitesting/random-graph-j-stats.rda")

rand.graphs = list()
for(i in 1:300){
  rand.graphs[[i]] <- randomGraph(maxnegtie = 0.1, pendchance = 0.02, badeggchance = 0.05)
}

rgk <- do.call('rbind', mclapply(rand.graphs, graphData, mc.cores=5))
save(rgk, file = "./piitesting/random-graph-k-stats.rda")


rand.graphs = list()
for(i in 1:300){
  rand.graphs[[i]] <- randomGraph(maxnegtie = 0.1, pendchance = 0.02, badeggchance = 0.1)
}

rgl <- do.call('rbind', mclapply(rand.graphs, graphData, mc.cores=5))
save(rgl, file = "./piitesting/random-graph-l-stats.rda")


table <- c()
for(i in 1:vcount(g)){
  for(j in 1:vcount(g)){
    if(i != j){
      xx <- polyroot((vec$pos[,i] - vec$neg[,i]) - (vec$pos[,j] - vec$neg[,j]))
      table <- c(table, Re(xx[round(Im(xx), digits = 12) == 0]))
    }
  }
}
table <- table[table < -0.5 & table > -0.9]
