negtie <- seq(0, 30, 10)
pendchance <- seq(0, .2, .05)
badegg <- seq(0, .1, .05)
combos <- expand.grid(negtie, pendchance, badegg) %>%
  rename(maxnegtie = Var1, pendchance = Var2, badeggchance = Var3)

combos.list <- split(combos, seq(nrow(combos)))
combos.list <- lapply(combos.list, as.list)

interruptable_saving <- function(num, point){
  for(n in num:60){
    interp_num <<- n
    for(i in point:100){
      interp_point <<- i
      print(paste(n, " : ", i, sep = ""))
      g <- do.call(randomGraph, combos.list[[n]])
      graph_data <- graphData(g)
      graph_data <-  cbind(graph_data, as.data.table(combos.list[[n]]))
      graph_data$step <- n
      graph_data$num <- i
      filename <- paste("./piitesting/random_graph_stats/", "step", n, "num", i, ".rda",sep = "_")
      save(g, graph_data, file = filename)
    }
    if(point!=1){
      point <- 1
    }
  }
}

interp_num = 1
interp_point = 1
interruptable_saving(interp_num, interp_point)

all.f <- list.files('./piitesting/random_graph_stats/')
data <- paste("./piitesting/random_graph_stats/", all.f[1], sep = "")
load(data)
all_graph_data <- graph_data
for(i in 2:length(all.f)){
  data <- paste("./piitesting/random_graph_stats/", all.f[i], sep = "")
  load(data)
  all_graph_data <- rbind(all_graph_data, graph_data)
}

no_neg_high <- all_graph_data[all_graph_data$epsStability == -0.98]
no_neg_low <- all_graph_data[all_graph_data$epsStability > -0.03] #==-0.02 giving false for values of -0.02
no_neg_mid <- all_graph_data[all_graph_data$epsStability == -0.5]
