#' Segmentation of left and right lungs from a CT scan
#'
#' @param img Original CT scan
#' @param lthresh Threshold used to find the lung and airways,
#' from the original CT scan. Anything below this threshold will
#' be used for the initial scan. Default: -300 HU.
#' @param verbose Print diagnostic messages.
#'
#' @return Lung mask, with left and right designation, of original CT scan.
#' Right lung has values = 1, left lung has values = 2, left/right lung overlap
#' (i.e. couldn't distinguish the left and right lungs) has values = 3, non-lung = 0.
#' @importFrom ANTsR maskImage
#' @importFrom ANTsRCore iMath labelClusters
#' @export
segment_lung_lr = function(img, lthresh = -300, verbose = TRUE){

  # Make all values positive, so 0s are 0s
  img = img + 1025
  lthresh = lthresh + 1025

  # Segmenting Lung and Airways
  lung_air_mask = segment_lung_airway(img, lthresh = lthresh, verbose = verbose)


  # Segmenting Airways
  air_mask = segment_airway(img, lung_air_mask, verbose)


  if (verbose) {
    message("# Removing Airways from Lung")
  }
  lung_mask = maskImage(lung_air_mask, 1-air_mask)
  img_masked = maskImage(img, lung_mask)
  img_mask = img_masked < 325
  lung_mask[img_mask == 0] = 0


  if (verbose) {
    message("# Finding Left and Right Lungs")
  }
  left_right_mask = labelClusters(lung_mask, minClusterSize = 50000)
  n_clus = length(unique(left_right_mask))

  # Left and right lung are still connected
  i = 1
  while(n_clus == 2){
    lung_mask = iMath(lung_mask, "ME", i)
    left_right_mask = labelClusters(lung_mask, minClusterSize = 50000)
    n_clus = length(unique(left_right_mask))
    mask1 = left_right_mask == 1
    mask1 = iMath(mask1, "MD", i)
    mask2 = left_right_mask == 2
    mask2 = iMath(mask2, "MD", i)
    left_right_mask = mask1 + mask2 * 2
    if(i > 5){break}
    i = i + 1
  }

  coord = sapply(1:n_clus, function(i){
    img = left_right_mask == i
    loc = which(as.array(img) == 1, arr.ind = T)
    x = median(loc[,1])
    return(x)
  })
  if (coord[1] < coord[2]){
    left_mask = left_right_mask == 1
    right_mask = left_right_mask == 2
  } else {
    left_mask = left_right_mask == 2
    right_mask = left_right_mask == 1
  }


  if (verbose) {
    message("# Smoothing")
  }
  left_mask = iMath(left_mask, "MD", 2)
  right_mask = iMath(right_mask, "MD", 2)
  left_right_mask = left_mask + right_mask * 2
  left_right_mask[left_right_mask == 3] = 0
  left_right_mask[lung_air_mask == 0] = 0

  return(left_right_mask)
}
