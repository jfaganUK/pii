#' Optimal Rank Beta
#'
#' This finds recommends a beta value between two extremes where the rank correlation between the two extremes best agrees.
#'
#' @param g An igraph graph object
#' @param comp.left The beta on the left for comparison
#' @param comp.right The beta on the right for comparison
#' @param starting.beta The beta value to start when starting the optimization
#' @examples
#' pii(g, pii.beta = optimal.rank.beta(g))
#'

optimal.rank.beta <-  function(g, comp.left = -0.9, comp.right = -0.5, starting.beta = -0.8) {
  rc.diff <- function(b, g, comp.left = comp.left, comp.right = comp.right) {
      p <- pii(g, pii.beta = b)
      pii.left <- pii(g, pii.beta = comp.left)
      pii.right <- pii(g, pii.beta = comp.right)
      rc.left = cor(p, pii.left, method = "spearman")
      rc.right = cor(p, pii.right, method = "spearman")
      return(abs(rc.left - rc.right))
    }

  o <- optim( par = c(starting.beta), fn = rc.diff, gr = NULL,
              g = g, comp.right = comp.right, comp.left = comp.left, method = 'Brent',
              lower = -1, upper = -0.001)
  cross.point <- o$par
  attr(cross.point, 'spearman.correlation') <- cross.y
  attr(cross.point, 'approx') <- o$value
  return(cross.point)
}
