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
      filename <- paste("./piitesting/random_graph_stats/", "step", n, "num", i, ".rda",sep = "_")
      save(graph_data, file = filename)
    }
    if(point!=1){
      point <- 1
    }
  }
}

interp_num = 1
interp_point = 1
interruptable_saving(interp_num, interp_point)
