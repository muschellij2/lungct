#' @title Coarse Masking of Image from Z-slicding
#' @description Calculates a Hull at each slice then stacks them
#'
#' @param x Binary \code{antsImage}
#' @param keep_slices number of slices to keep if needed
#'
#' @return Array or \code{antsImage} object
#' @export
#' @importFrom ANTsRCore is.antsImage
#' @importFrom ptinpoly pip2d
coarse_body = function(x, keep_slices = 20) {
  
  bb = as.array(x) > 0
  dimg = dim(bb)
  number_slices = dimg[3]
  dslice = dimg[1:2]
  ind = which(bb, arr.ind = TRUE)
  omat = array(FALSE, dim = dimg)
  
  
  ind = as.data.frame(ind)
  ps = split(ind, ind[,3])
  ps = lapply(ps, function(x){
    as.matrix(x[, 1:2])
  })
  
  # get values that are true for ellipse
  # same for all slices
  L = lapply(dslice, seq)
  test.ind = as.matrix(expand.grid(L))
  
  res = lapply(ps, ptinpoly::pip2d, Queries = test.ind)
  
  i = 1
  for (i in seq(number_slices)) {
    # print(i)
    r = res[[i]] >= 0
    get_ind = test.ind[r, ]
    omat[,,i][get_ind] = TRUE
  }
  
  prob = apply(omat, c(1,2), mean)
  prob = prob > (keep_slices/number_slices)
  
  for (i in seq(number_slices)) {
    # print(i)
    omat[,,i][!prob] = FALSE
  }
  if (is.antsImage(x)) {
    omat = as.antsImage(omat, reference = x)
  }
  return(omat)
}
