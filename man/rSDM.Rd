% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparseDirichletMixture.R
\name{rSDM}
\alias{rSDM}
\title{Random draw for probability vector under sparse Dirichlet mixture model}
\usage{
rSDM(alpha, beta, logscale = FALSE, allout = FALSE)
}
\description{
Random draw for probability vector under sparse Dirichlet mixture model
}
\examples{
K = 5 # number of categories
a = rep(1.0, K) # Dirichlet shape parameter vector
b = 10.0 # SDM parameter
rSDM(a, b) # single draw from prior
x = rbinom(K, 10, 0.3) # random count data
rSDM(a+x, b) # single draw from posterior
}
