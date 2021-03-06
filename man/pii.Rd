% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pii.r
\name{pii}
\alias{pii}
\title{Political Independence Index}
\usage{
pii(g, pii.beta = -0.8, e.dist = NULL, triadic = F, pii.delta = 0.1,
  max.degree = NULL)
}
\arguments{
\item{g}{An igraph graph}

\item{pii.beta}{Should the vertex and edge names be added to the rows and columns of the matrix.}

\item{e.dist}{(optional) an edge.distance matrix, if calculated ahead of time.}
}
\description{
Returns a numeric vector of PII scores for each vertex in network.
}
\examples{
pii(g)
}
