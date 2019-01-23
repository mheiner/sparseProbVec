#' Convert Dirichlet shape parameter vector into Generalized Dirichlet beta parameters
#' @export
shape_Dir2genDir = function(a) {
  K = length(a)
  rcra = rev( cumsum( rev(a[-1]) ) )
  cbind( a=a[-K], b=rcra )
}

#' SBM correction to Generalized Dirichlet parameters
#' @export
delCorrection_SBM = function(p_small, p_large, K) {
  if (p_large == 0.0) {
    out = ((1.0 - p_small)*(K - 1) + 1) / K
  } else {
    aa = 1.0 - p_large
    bb = 1.0 - aa^K
    out = ((1.0 - p_small)*bb/p_large + p_small) / K
  }
  out
}

#' Simulate from sparse stick-breaking mixture prior
#' @export
#' @examples
#' a = c(5, 2, 3)
#' K = length(a)
#' rSBM(K=K, p_small=0.5, p_large=0.0, eta=1.0e3, gam=1.5, del=1.5)
#' gamdel = shape_Dir2genDir(a) # tapered parameters to mimic Dirichlet special case
#' rSBM(K=K, p_small=0.5, p_large=0.0, eta=1.0e3, gam=gamdel[,1], del=gamdel[,2]) # tapered to mimic Dirichlet special case
rSBM = function(K, p_small, p_large=0.0, eta, gam, del, logout=FALSE, logcompute=TRUE) {
  if (!logcompute && logout) stop("logout=TRUE requires logcompute=TRUE")

  stopifnot(eta > 1.0)
  n = K - 1

  if ( length(p_small) == 1) {
    p_small = rep(p_small, n)
  }

  if ( length(p_large) == 1) {
    p_large = rep(p_large, n)
  }

  ## proportions
  p_GenDir = 1.0 - p_small - p_large
  stopifnot(all(p_small >= 0.0), all(p_large >= 0.0), all(p_GenDir >= 0.0), all(p_GenDir <= 1.0))


  if ( length(gam) == 1 ) {
    gam = rep(gam, n)
  }

  if ( length(del) == 1 ) {
    del = rep(del, n)
  }


  ## group membership
  xi = numeric(n)
  for ( i in 1:n ) {
    xi[i] = sample.int(n=3, size=1, replace=TRUE, prob=c(p_small[i], p_GenDir[i], p_large[i]))
  }

  ## latent z
  z = numeric(n)
  lz = numeric(n)
  loneminusz = numeric(n)

  if ( logcompute ) {
    for ( i in 1:n ) {

      if ( xi[i] == 1 ) {
        lx1 = log( rexp(1, 1.0) )
        lx2 = log (rgamma(1, eta, 1.0) )
      } else if ( xi[i] == 2 ) {
        lx1 = log( rgamma(1, gam[i], 1.0) )
        lx2 = log( rgamma(1, del[i], 1.0) )
      } else {
        lx1 = log (rgamma(1, eta, 1.0) )
        lx2 = log( rexp(1, 1.0) )
      }

      lxm = max(lx1, lx2)
      lxsum = lxm + log( exp(lx1 - lxm) + exp(lx2 - lxm) )
      lz[i] = lx1 - lxsum
      loneminusz[i] = lx2 - lxsum

    }
  } else { # not logcompute
    for ( i in 1:n ) {
      if ( xi[i] == 1 ) {
        z[i] = rbeta(1, 1.0, eta)
      } else if ( xi[i] == 2 ) {
        z[i] = rbeta(1, gam[i], del[i])
      } else {
        z[i] = rbeta(1, eta, 1.0)
      }
    }
  }

  if (logcompute) {
    lw = numeric(K)
    lwhatsleft = 0.0
    for ( i in 1:n ) {
      lw[i] = lz[i] + lwhatsleft
      # lwhatsleft = lwhatsleft + log( 1.0 - exp(lw[i] - lwhatsleft) ) # logsumexp
      lwhatsleft = lwhatsleft + loneminusz[i]
    }
    lw[K] = lwhatsleft
    w = exp(lw)
  } else {
    w = numeric(K)
    whatsleft = 1.0
    for ( i in 1:n ) {
      w[i] = z[i] * whatsleft
      whatsleft = whatsleft - w[i]
    }
    w[K] = whatsleft
  }

  stopifnot(all(w > 0.0), all(w < 1.0))

  if (logout) {
    return( lw )
  } else {
    return( w )
  }
}




#' Sample Conjugte Posterior: Probability vector from multinomial counts under sparse stick-breaking mixture prior
#'
#' @param x real vector, counts in data
#' @return list
#' @export
#' @examples
#' x = c(5, 2, 3)
#' K = length(x)
#' rPostSBM(x=x, p_small=0.5, p_large=0.0, eta=1.0e3, gam=1.5, del=1.5)
rPostSBM = function(x, p_small, p_large=0.0, eta, gam, del, w_logout=FALSE, z_logout=FALSE) {

  stopifnot(eta > 1.0)

  ## Gibbs draw for stick breaking betas and resulting weights
  ## based on Generalized Dirichlet Dist.
  K = length(x)
  n = K - 1


  if ( length(p_small) == 1 ) {
    p_small = rep(p_small, n)
  }

  if ( length(p_large) == 1 ) {
    p_large = rep(p_large, n)
  }

  if ( length(del) == 1 ) {
    del = rep(del, n)
  }

  rcrx = rev( cumsum( rev(x[-1]) ) )

  ## mixture component 1 update
  a_1 = 1.0 + x[1:n]
  b_1 = eta + rcrx

  ## mixture component 2 update
  a_2 = gam + x[1:n]
  b_2 = del + rcrx

  ## mixture component 3 update
  a_3 = eta + x[1:n]
  b_3 = 1.0 + rcrx

  aa = list(a_1, a_2, a_3)
  bb = list(b_1, b_2, b_3)

  ## proportions
  p_GenDir = 1.0 - p_small - p_large
  stopifnot(all(p_small >= 0.0), all(p_large >= 0.0), all(p_GenDir >= 0.0), all(p_GenDir <= 1.0))


  ## calculate posterior mixture weights
  leta = log(eta)
  logweight1 = log(p_small) + leta + lbeta(a_1, b_1)
  logweight2 = log(p_GenDir) - lbeta(gam, del) + lbeta(a_2, b_2)
  logweight3 = log(p_large) + leta + lbeta(a_3, b_3)

  lmax = pmax(logweight1, logweight2, logweight3)
  ldenom = lmax + log( exp(logweight1 - lmax) +
                         exp(logweight2 - lmax) +
                         exp(logweight3 - lmax) ) # logsumexp

  lp_xi = cbind(logweight1 - ldenom, logweight2 - ldenom, logweight3 - ldenom)

  ## draw mixture membership indicators
  xi = numeric(n)
  for ( i in 1:n ) {
    xi[i] = sample.int(n=3, size=1, prob=exp(lp_xi[i,]))
  }

  grp_indx = list()
  grp_n = list()
  lz = numeric(n)
  loneminusz = numeric(n)

  for ( ii in 1:3 ) {
    grp_indx[[ii]] = which(xi == ii)
    grp_n[[ii]] = length(grp_indx[[ii]])
    if ( grp_n[[ii]] > 0 ) {
      lx1 = log( rgamma(grp_n[[ii]], aa[[ii]][grp_indx[[ii]]], 1.0) )
      lx2 = log( rgamma(grp_n[[ii]], bb[[ii]][grp_indx[[ii]]], 1.0) )
      lxm = pmax(lx1, lx2)
      lxsum = lxm + log( exp(lx1 - lxm) + exp(lx2 - lxm) )
      lz[grp_indx[[ii]]] = lx1 - lxsum
      loneminusz[grp_indx[[ii]]] = lx2 - lxsum
    }
  }

  ## replaced by the preceding loop
  # if (n1 > 0) {
  #   lx1 = log( rgamma(n1, a_1[grp1indx], 1.0) )
  #   lx2 = log( rgamma(n1, b_1[grp1indx], 1.0) )
  #   lxm = pmax(lx1, lx2)
  #   lxsum = lxm + log( exp(lx1 - lxm) + exp(lx2 - lxm) )
  #   lz[grp1indx] = lx1 - lxsum
  # }
  # if (n2 > 0) {
  #   lx1 = log( rgamma(n2, a_2[grp2indx], 1.0) )
  #   lx2 = log( rgamma(n2, b_2[grp2indx], 1.0) )
  #   lxm = pmax(lx1, lx2)
  #   lxsum = lxm + log( exp(lx1 - lxm) + exp(lx2 - lxm) ) # logsumexp
  #   lz[grp2indx] = lx1 - lxsum
  # }
  # if (n3 > 0) {
  #   lx1 = log( rgamma(n3, a_3[grp3indx], 1.0) )
  #   lx2 = log( rgamma(n3, b_3[grp3indx], 1.0) )
  #   lxm = pmax(lx1, lx2)
  #   lxsum = lxm + log( exp(lx1 - lxm) + exp(lx2 - lxm) ) # logsumexp
  #   lz[grp3indx] = lx1 - lxsum
  # }

  ## break the stick

  ## not on log scale
  # w = numeric(K)
  # whatsleft = 1.0
  # for ( i in 1:n ) {
  #   w[i] = z[i] * whatsleft
  #   whatsleft = whatsleft - w[i]
  # }
  # w[K] = whatsleft

  ## on log scale
  lw = numeric(K)
  lwhatsleft = 0.0
  for ( i in 1:n ) {
    lw[i] = lz[i] + lwhatsleft
    # lwhatsleft = lwhatsleft + log( 1.0 - exp(lw[i] - lwhatsleft) ) # logsumexp
    lwhatsleft = lwhatsleft + loneminusz[i]
  }
  lw[K] = lwhatsleft

  if (w_logout) {
    out_w = lw
  } else {
    out_w = exp(lw)
  }

  if (z_logout) {
    out_z = lz
  } else {
    out_z = exp(lz)
  }

  list(w=out_w, z=out_z, xi=xi)
}
