#' @title Coarse Masking of Image from Z-slicding
#' @description Calculates an Ellipsoid at each slice then stacks them
#'
#' @param x Binary \code{antsImage}
#'
#' @return Array or \code{antsImage} object
#' @export
#' @importFrom ANTsR is.antsImage
#' @importFrom stats cov
ellipsoid_body = function(x) {
  bb = as.array(x)
  dimg = dim(bb)
  number_slices = dimg[3]
  dslice = dimg[1:2]
  omat = array(FALSE, dim = dimg)

  # get values that are true for ellipse
  ind = which(bb > 0, arr.ind = TRUE)
  L = lapply(dslice, seq)
  # same for all slices
  X = as.matrix(expand.grid(L))

  # need df for split
  ind = as.data.frame(ind)
  ps = split(ind, ind[,3])
  ps = lapply(ps, function(x){
    as.matrix(x[, 1:2])
  })
  # Getting mu/sigma for ellipse
  vals = lapply(ps, function(x){
    L = list(mu = colMeans(x),
             Sigma = cov(x))
  })
  # Taken from SIBER, but faster because no looping over solve
  pointsToEllipsoid2 = function(X, Sigma, mu) {
    eig <- eigen(Sigma)
    SigSqrt = eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
    xtx = t(SigSqrt) %*% SigSqrt
    ss = solve(xtx)
    xx = t(t(X) - mu)
    Z = ss %*% t(SigSqrt) %*% t(xx)
    Z = t(Z)
  }
  # Taken from SIBER
  ellipseInOut2 = function(Z, p = 0.95, radius = NULL) {
    if (is.null(radius)) {
      radius <- stats::qchisq(p, df = ncol(Z))
    }
    inside <- rowSums(Z^2) < radius
    return(inside)
  }
  for ( i in seq(number_slices)) {
    x = vals[[i]]
    Z = pointsToEllipsoid2(X, Sigma = x$Sigma, mu = x$mu)
    inside = ellipseInOut2(Z = Z, p = 0.95)
    mat = array(inside, dim = dslice)
    omat[,,i] = mat
  }
  if (is.antsImage(x)) {
    omat = as.antsImage(omat, reference = x)
  }
  return(omat)
}
