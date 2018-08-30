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

  orig_img = check_ants(img)
  img = antsImageClone(img)
  if (verbose) {
    message("# Resampling Image to 1x1x1")
  }
  img = resampleImage(img, c(1,1,1))
  # Make all values positive, so 0s are 0s
  #img = img + 1025
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

  if (verbose) {
    message("# Finding Left and Right Lungs")
  }
  left_right_mask = labelClusters(lung_mask, minClusterSize = 50000)
  n_clus = length(unique(left_right_mask))

  if (n_clus == 1) {
    message("# Error: Can't find lungs, returning current mask")
    return(lung_mask)
  }

  # Left and right lung are still connected
  i = 0
  lung_mask2 = antsImageClone(lung_mask)
  while(n_clus == 2){
    i = i + 1
    new_lthresh = lthresh - 25*i
    if(new_lthresh < 0){
      message("# Error: Can't distinguish left/right lungs, returning left/right combined mask")
      return(lung_mask)
      }
    img_mask = img_masked < new_lthresh
    lung_mask2[img_mask == 0] = 0
    left_right_mask = labelClusters(lung_mask2, minClusterSize = 50000)
    n_clus = length(unique(left_right_mask))
  }

  # Correct number of clusters
  if (n_clus == 3) {

    # Find coordinates of left and right clusters
    coord = sapply(1:n_clus, function(i){
      img = left_right_mask == i
      loc = which(as.array(img) == 1, arr.ind = T)
      x = median(loc[,1])
      return(x)
    })
    if (coord[1] < coord[2]){
      right_mask = left_right_mask == 1
      left_mask = left_right_mask == 2
    } else {
      right_mask = left_right_mask == 2
      left_mask = left_right_mask == 1
    }


    if (verbose) {
      message("# Finishing Touches")
    }
    left_mask = iMath(left_mask, "MC", 3)
    right_mask = iMath(right_mask, "MC", 3)
    left_right_mask = right_mask + left_mask * 2
    left_right_mask[left_right_mask == 3] = 0
    left_right_mask[lung_air_mask == 0] = 0

  } else{message("# Error: Too many clusters")}


  if (verbose) {
    message("# Resampling Back to Original Image Size")
  }
  left_right_mask = resampleImage(left_right_mask,
                            resampleParams = dim(orig_img),
                            useVoxels = TRUE,
                            interpType = 1)
  return(left_right_mask)

}

























