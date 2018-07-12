#' @rdname segment_lung
#' @export
#' @importFrom ANTsR connectedThreshold maskImage
#' @importFrom ANTsRCore iMath
quick_lung_mask = function(img,
                           lthresh = -300) {
  img = check_ants(img)

  # just get rid of the standard background just in case
  bg_mask = iMath(img < -950, "GetLargestComponent")
  mask = maskImage(img > -1023 & img < -900, 1 - bg_mask)
  ind = which(as.array(mask) > 0, arr.ind = TRUE)
  seed = floor(colMeans(ind))
  d = colMeans((t(ind) - seed)^2)
  seed = ind[which.min(d),]

  simg = connectedThreshold(img, seed = seed, upper = lthresh, lower = -1023)

}