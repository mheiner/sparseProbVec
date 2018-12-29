
#' Log Sum Exp function
#'
#' Numerically stable computation of log(sum(exp(x))).
#'
# ' @keywords
# ' @seealso \code{\link{diag_vech_indx}}
#' @param x real vector, the arguments to exp(). If log(sum(y)) is desired, supply x=log(y).
#' @return real number, numerically stable log(sum(exp(x)))
#' @export
#' @examples
#' logsumexp(1:5)
#' log(sum(exp(1:5)))
logsumexp <- function(x) {
  m = max(x)
  m + log(sum(exp(x-m)))
}



#' Single draw from Dirichlet distribution
#'
#' @keywords Dirichlet
# ' @seealso \code{\link{}}
#' @param alpha real vector, shape parameter.
#' @param logscale, should probability vector be returned on the log scale?
#' @return real probability vector, single draw from the Dirichlet(alpha) distribution
#' @export
#' @examples
#' rDiriclet(rep(1,5))
rDirichlet = function(alpha, logscale=FALSE) {
  len = length(alpha)
  xx = rgamma(len, alpha, 1.0)
  if (logscale) {
    lx = log(xx)
    out = lx - sparseProbVec::logsumexp(lx)
  } else {
    out = xx / sum(xx)
  }
  out
}
