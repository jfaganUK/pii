### Mining the simulation results ##############################################
rm(list=ls())
gc()

### Libraries
library(glmnet)
library(data.table)
library(dplyr)
library(ggplot2)
library(texreg)

### Load the data
load('all-random-graph-stats.rda')
w <- allringgraphs; rm('allringgraphs')
w <- na.omit(w)

### The two variable sets
xvars <- c("meanDegree","density","diameter","nodeCount","edgeCount","degCentralization","avgPathLength","numPosEdge","numNegEdge","propNegEdge","modularity","meanTrans","avgMinDistToNegEdge","avgDistOfNegEdge")
yvars <- c('rankCor', 'lowCorBeta')

fixVarName <- function(ss) {
  name.change <- data.table(xvars = xvars %>% sort,
                            new.var = c('Avg. Distance to Neg. Edge', 'Avg. Min. Dist. to Neg. Edge', 'Avg. Path Length',
                                        'Degree Centralization', 'Density', 'Density', 'Diameter',
                                        'Edge Count', 'Mean Degree', 'Mean Transitivity', 'Modularity',
                                        'Node Count', 'Number of Neg. Edge', 'Number of Pos. Edge', 'Proportion of Neg. Edges'))
  for(i in 1:nrow(name.change)) {
    ss <- gsub(name.change$xvars[i], name.change$new.var[i], ss)
  }
  ss <- gsub('X', ' * ', ss)
  return(ss)
}

corstarsl <- function(x){
  require(Hmisc)
  x <- as.matrix(x)
  R <- rcorr(x)$r
  p <- rcorr(x)$P

  ## define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))

  ## trunctuate the matrix that holds the correlations to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]

  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")

  ## remove upper triangle
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew)

  ## remove last column and return the matrix (which is now a data frame)
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  return(Rnew)
}



### Chart theme
thm <- theme(
  panel.background = element_rect(fill = 'white'),
  plot.background = element_rect(fill = 'white'),
  legend.background = element_rect(fill = 'white'),
  axis.ticks = element_blank(),
  axis.text = element_text(family = 'Futura Medium'),
  axis.title = element_text(family = 'Futura Medium', size = 15),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(color = '#AAAAAA')
)

### A plot of, what I think, is the primary variable.
ggplot(w[avgMinDistToNegEdge > 0], aes(y=rankCor, x=avgMinDistToNegEdge)) +
  scale_x_continuous('Avg. Distance to First Negative Edge') +
  scale_y_continuous('Rank Correlation\n(PII of beta = 0.5 & beta = 0.9)') +
  geom_point() +
  geom_smooth(method='lm') +
  thm

### Stats table ##########################################################

ss <- w[, c('rankCor', xvars), with = F] %>%
  melt %>% group_by(variable) %>%
  summarize(mean = mean(value), sd = sd(value), min = min(value), max = max(value)) %>%
  mutate(`Mean (SD) [Min, Max]` = sprintf('%3.2f (%1.2f) [%1.2f, %1.2f]', mean, sd, min, max)) %>%
  select(variable, `Mean (SD) [Min, Max]`)
library(formattable)
formattable(ss %>% data.frame)

corstarsl(w[, c('rankCor', xvars), with = F] ) %>%
  write.csv('feature_correlation_table.csv')

### Scaling and Creating Terms #################################################
## note, I don't want to scale shit right now
for(cc in xvars)
  w[, eval(cc) := scale(.SD), .SDcols = cc]


f <- glmnet(w %>% select(match(xvars, names(w))) %>% as.matrix,
            w$rankCor, alpha = 1)
coef(f, s = 0.005)

### Create interaction terms
interaction.list <- combn(xvars, 2) %>% t %>% data.frame
interaction.list <- rbind(interaction.list, data.frame(X1 = xvars, X2 = xvars))
createInteractionVariable <- function(il) {
  v1 <- il['X1']
  v2 <- il['X2']
  w1 <- w[, v1, with=F]
  w2 <- w[, v2, with=F]
  w[, eval(paste0(v1,'X',v2)) := w1*w2]
}
invisible({
  apply(interaction.list, 1, createInteractionVariable)
})

### Build datasets
y <- as.matrix(w[, rankCor])
nm <- names(w)
nm <- nm[nm != 'graphid']
nm <- nm[nm != 'rankCor']
nm <- nm[nm != 'lowCorBeta']
x <- as.matrix(w[,nm, with=F])
colnames(x) <- nm

### Penalized Regression #######################################################
### Fit the model
fit1 <- cv.glmnet(x, y, alpha=0.2)  # elasticnet
fit2 <- cv.glmnet(x, y, alpha=0)    # ridge regression
fit3 <- cv.glmnet(x, y, alpha=1)    # LASSO

### Get the results
f <- fit1
plot(f)
coef(f, s='lambda.1se') %>%
  as.matrix %>% data.table(keep.rownames=T) %>% rename(`coef` = `1`) %>%
  filter(abs(coef) > 0) %>% nrow

which(f$nzero == 13)[1]

### get the coef table
fit.coef <-  coef(f, s='lambda.1se') %>%
  as.matrix %>% data.table(keep.rownames=T) %>% rename(`coef` = `1`) %>%
  filter(abs(coef) > 5e-3) %>%  filter(rn != '(Intercept)') %>%
  arrange(coef) %>% mutate(rn = factor(rn, levels=.$rn))

### what proportion of the total beta is accounted for
sum(abs(fit.coef$coef)) / (abs(coef(f, s='lambda.1se')[,1]) %>% sum)

### plot
ggplot(fit.coef, aes(y=rn, x=coef)) +
  geom_point() +
  thm

### table
library(formattable)
o <- fit.coef %>% rename(`Variable` = rn, Beta = coef) %>%
  mutate(Beta = round(Beta, 3))
formattable(o)

### Count the number of times variable shows up in a coefficient
rn <- coef(f, s='lambda.1se') %>%
  as.matrix %>% data.table(keep.rownames=T) %>% rename(`coef` = `1`) %>%
  filter(abs(coef) > 0) %>% .[['rn']]
xvar.imp <- numeric(length(xvars))
names(xvar.imp) <- xvars
for(i in xvars)
{
  xvar.imp[i] <- grepl(i, rn) %>% sum
}

fit.coef

### Random Forest ##############################################################
library(randomForest)
set.seed(10101982)
f.rf1 <- randomForest(rankCor ~ meanDegree + density + diameter + nodeCount + edgeCount + degCentralization + avgPathLength + numPosEdge + numNegEdge + propNegEdge + modularity + meanTrans + avgMinDistToNegEdge + avgDistOfNegEdge,
                      data = w, importance = T, prox = T)
cor(w$rankCor, predict(f.rf1))^2
plot(w$rankCor, predict(f.rf1))
o <- round(importance(f.rf1), 2) %>% data.table(keep.rownames=T) %>%
  arrange(desc(`%IncMSE`))
library(formattable)
formattable(o)

nd.density <- data.table(meanDegree = 0, density = 0,
                         diameter = 0, nodeCount = 0, edgeCount = 0, degCentralization = 0,
                         avgPathLength = 0, numPosEdge = 0, numNegEdge = 0,
                         propNegEdge = c(rep(min(w$propNegEdge), 100), rep(mean(w$propNegEdge), 100), rep(max(w$propNegEdge), 100)),
                         modularity = 0, meanTrans = 0,
                         avgMinDistToNegEdge = rep(seq(min(w$avgMinDistToNegEdge), max(w$avgMinDistToNegEdge), length.out = 100), 3),
                         avgDistOfNegEdge = 0)

p <- predict(f.rf1, newdata =  nd.density)
d <- cbind(nd.density, p)

ggplot(d, aes(avgMinDistToNegEdge, p, group = propNegEdge, color = factor(propNegEdge, labels = c('Low', 'Med', 'High')))) + geom_line() +
  scale_color_discrete('Prop. Neg. Edges') +
  scale_x_continuous('Avg. Distance to First Negative Edge') +
  scale_y_continuous('Spearman Correlation\n(Random Forest Est.)') +
  thm

d <- importance(f.rf1) %>% data.table(keep.rownames=T) %>% mutate(rn = fixVarName(rn))
ggplot(d, aes(x = `%IncMSE`, y = `IncNodePurity`, label = rn)) +
  scale_x_continuous('% Increase in MSE') +
  scale_y_continuous('Increase in Node Purity') +
  geom_text() + thm


library(tsne)
f.rf1.tsne <- tsne(f.rf1$proximity)



### OLS ########################################################################
# lm1 <- lm(rankCor ~ meanDegree + density + diameter + nodeCount + edgeCount + degCentralization + avgPathLength + numPosEdge + numNegEdge + propNegEdge + modularity + meanTrans + avgMinDistToNegEdge + avgDistOfNegEdge, data = w)
lm1 <- lm(rankCor ~ propNegEdge + avgMinDistToNegEdge + avgPathLength, data = w)
summary(lm1)
screenreg(lm1)

library(car)
vif(lm1)

pp <- data.table(original = y, ols = predict(lm1), elastic = predict(fit1, newx = x), rand.forest = predict(f.rf1))
cor(pp)^2

lm2 <- lm(rankCor ~ avgMinDistToNegEdge + propNegEdge  + meanDegree + nodeCount + meanDegree:nodeCount +
            avgMinDistToNegEdgeXavgMinDistToNegEdge + avgMinDistToNegEdgeXpropNegEdge, data=w)
screenreg(lm1, custom.coef.names = c('Intercept', 'Avg. Dist. to First Neg. Edge', 'Proportion of Neg. Edges', 'Mean Degree', 'Size', 'ADFNE^2', 'ADFNE * PropNegEdge', 'MeanDegree * Size'))
dev.new()
par(bg = '#eeeeee')
plotreg(lm1, insignif.light = '#999999',insignif.medium='#777777',insignif.dark='#555555',
        signif.light = "#c5dbe9", signif.medium = "#5a9ecc", signif.dark = "#1c5ba6",
        xlim =c(-0.1,0.1), custom.model.names = 'Beta Rank Correlation',
        custom.coef.names = c('Intercept', 'Avg. Dist. to First Neg. Edge', 'Proportion of Neg. Edges', 'Mean Degree', 'Size', 'ADFNE^2', 'ADFNE * PropNegEdge', 'MeanDegree * Size'))
dev.off()

cv <- data.table(rc1 = w$rankCor, rc2 = predict(fit1, x))
cor(cv$rc1, cv$rc2)
ggplot(cv, aes(rc2.1, rc1)) + geom_point() +
  scale_y_continuous('Empirical Rank Correlation') +
  scale_x_continuous('Predicted Rank Correlation') +
  geom_smooth(method='lm') + thm

