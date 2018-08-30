#' Segmentation of lung and airway from a CT scan
#'
#' @param img Original CT scan
#' @param lthresh Threshold used to find the lung and airways,
#' from the original CT scan. Anything below this threshold will
#' be used for the initial scan. Default: -300 HU.
#' @param verbose Print diagnostic messages.
#'
#' @return Lung mask of lung and airways from original CT scan.
#' Mask has values = 1, background = 0.
#' @importFrom ANTsR maskImage
#' @importFrom ANTsRCore iMath labelClusters antsImageClone
#' @export
segment_lung_airway = function(
  img, lthresh = -300, verbose = TRUE){

  # Simple thresholding
  if (verbose) {
    message("# Segmenting Lung and Airways: Thresholding")
  }
  body_img = img <= lthresh


  # Find the lungs and trachea from the image
  if (verbose) {
    message("# Segmenting Lung and Airways: Largest Component")
  }
  first_img = iMath(img = body_img, operation = "GetLargestComponent")

  # Check to see if we grabbed the background instead of the lung/trachea
  # The back 10 slice and/or front 10 slices should have mean ~ 1, if selected
  n_max = dim(first_img)[2]
  mean_back = mean(apply(as.array(first_img)[,1:10,],2,mean))
  mean_front = mean(apply(as.array(first_img)[,(n_max-9):n_max,],2,mean))

  # Remove background, as necessary
  if(mean_back > .5 | mean_front > .5)
  {
    first_img_hold = antsImageClone(first_img)

    # Remove background from rest of image
    first_img = (1-first_img_hold) * body_img
    first_img = iMath(img = first_img, operation = "GetLargestComponent")

    # Check to see if that fixed the problem
    mean_back = mean(apply(as.array(first_img)[,1:10,],2,mean))
    mean_front = mean(apply(as.array(first_img)[,(n_max-9):n_max,],2,mean))

    # Keep fixing, as necessary
    if(mean_back > .5 | mean_front > .5) # Still background
    {
      first_img_hold[first_img==1] = 1
      first_img = (1-first_img_hold) * body_img
      first_img = iMath(img = first_img, operation = "GetLargestComponent")

      # Check to see if that fixed the problem
      mean_back = mean(apply(as.array(first_img)[,1:10,],2,mean))
      mean_front = mean(apply(as.array(first_img)[,(n_max-9):n_max,],2,mean))
      if(mean_back > .5 | mean_front > .5){stop("Can't find lungs beneath background noise. More coding needed")}

    }
  }
  body_img = iMath(first_img, "MC", 2)

  return(body_img)
}
