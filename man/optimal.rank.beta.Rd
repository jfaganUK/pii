% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimal.rank.beta.r
\name{optimal.rank.beta}
\alias{optimal.rank.beta}
\title{Optimal Rank Beta}
\usage{
optimal.rank.beta(g, comp.left = -0.9, comp.right = -0.5,
  starting.beta = -0.8)
}
\arguments{
\item{g}{An igraph graph object}

\item{comp.left}{The beta on the left for comparison}

\item{comp.right}{The beta on the right for comparison}

\item{starting.beta}{The beta value to start when starting the optimization}
}
\description{
This finds recommends a beta value between two extremes where the rank correlation between the two extremes best agrees.
}
\examples{
pii(g, pii.beta = optimal.rank.beta(g))

}
