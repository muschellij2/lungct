#' Segmentation of airways from a CT scan
#'
#' @param img Original CT scan
#' @param lung_air_mask Mask of original CT scan from segment_lung_airway
#' @param verbose Print diagnostic messages
#'
#' @return Mask of airways only
#' @export
#' @importFrom ANTsR maskImage
#' @importFrom ANTsRCore antsImageClone labelClusters iMath
#' @importFrom stats quantile
segment_airway = function(img, lung_air_mask, verbose = TRUE) {
  if (verbose) {
    message("# Segmenting Airways: Thresholding")
  }
  lung_air_mask = check_ants(lung_air_mask)
  img = check_ants(img)

  air_mask = antsImageClone(lung_air_mask)
  air = maskImage(img,lung_air_mask)
  air = air[air>0]
  air_thres = quantile(air,.05)
  air_mask[img > air_thres] = 0
  air_mask = iMath(air_mask, "MO", 1)
  air_mask = iMath(air_mask, "MD", 1)

  if (verbose) {
    message("# Segmenting Airways: Connected Components")
  }
  air_mask = labelClusters(air_mask, minClusterSize = 10000)
  n_clus = length(unique(air_mask))
  if(n_clus == 1){
    air_mask = antsImageClone(lung_air_mask)
    air_mask[img > air_thres*.5] = 0
    air_mask = iMath(air_mask, "MO", 1)
    air_mask = iMath(air_mask, "MD", 1)
    air_mask = labelClusters(air_mask, minClusterSize = 10000)
    n_clus = length(unique(air_mask))
  }
  if(n_clus > 2){
    coord = sapply(1:n_clus, function(i){
      iimg = air_mask == i
      loc = which(as.array(iimg) == 1, arr.ind = T)
      x = median(loc[,1])
      return(x)
    })
    center = dim(air_mask)[1]/2
    loss = abs(center - coord)
    clus = which(loss == min(loss, na.rm = T))
    air_mask = air_mask == clus
  }

  if (verbose) {
    message("# Segmenting Airways: Smoothing")
  }
  air_mask = iMath(air_mask, "MD", 1)
  return(air_mask)
}
