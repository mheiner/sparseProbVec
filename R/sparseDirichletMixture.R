#' Random draw for probability vector under sparse Dirichlet mixture model
#'
#'
# ' @keywords conjugate Gibbs normal
# ' @seealso \code{\link{}}
# ' @param x real vector, counts in data
# ' @param alpha real (numeric or vector \code{length(x) - 1}), defaults to 1.0 (approximating DP stick weights)
# ' @param beta real, (numcric or vector \code{length(x) - 1}), the alpha in approximating DP stick weights
# ' @param v0 real, prior variance
# ' @param nsim int, number of samples to draw, defaults to 1
# ' @return real, draw from the normal posterior
#' @export
#' @examples
#' K = 5 # number of categories
#' a = rep(1.0, K) # Dirichlet shape parameter vector
#' b = 10.0 # SDM parameter
#' rSDM(a, b) # single draw from prior
#' x = rbinom(K, 10, 0.3) # random count data
#' rSDM(a+x, b) # single draw from posterior
rSDM = function(alpha, beta, logscale=FALSE, allout=FALSE) {
  K = length(alpha)
  X = matrix(rep(alpha, each=K), nrow=K) + beta*diag(K)
  lgMat = lgamma(X)
  lpg = rowSums(lgMat)
  lpg_denom = logsumexp(lpg) # not even necessary, could just exponentiate lpgs and use weights
  w = exp(lpg - lpg_denom)

  z = sample.int(K, size=1, prob=w)
  theta_out = sparseProbVec::rDirichlet(X[z,], logscale)

  if ( allout ) {
    return( list(theta=theta_out, z=z) )
  } else {
    return(theta_out)
  }
}
