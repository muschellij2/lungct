#' Segmentation of left and right lungs from a CT scan
#'
#' @param img Original CT scan
#' @param lthresh Threshold used to find the lung and airways,
#' from the original CT scan. Anything below this threshold will
#' be used for the initial scan. Default: -300 HU.
#' @param verbose Print diagnostic messages.
#'
#' @return Lung mask, with left and right designation, of original CT scan.
#' Right lung has values = 1, left lung has values = 2, non-lung = 0.
#' @importFrom ANTsR maskImage
#' @importFrom ANTsRCore iMath labelClusters makeImage
#' @export
segment_lung_lr = function(img, lthresh = -300, verbose = TRUE){

    orig_img = check_ants(img)
    img = antsImageClone(orig_img)
    if (verbose) {
      message("# Resampling Image to 1x1x1")
    }
    img = resampleImage(img, c(1,1,1))


    # Make all values positive, so 0s are 0s
    img = img + 1025
    orig_lthresh = lthresh
    lthresh = lthresh + 1025


    # Segmenting Lung and Airways (using code from segment_lung)
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
        if(mean_back > .5 | mean_front > .5){
          first_img = segment_lung(orig_img, lthresh = orig_lthresh)
          stop("Can't find lungs beneath background noise. More coding needed")
          }

      }
    }
    lung_air_mask = iMath(first_img, "MC", 2)



    # Segmenting Airways
    if (verbose) {
      message("# Segmenting Airways: Thresholding")
    }
    air_mask = antsImageClone(lung_air_mask)
    air = maskImage(img,lung_air_mask)
    air_thres = quantile(air[air>0],.04)
    air_mask[img > air_thres] = 0

    if (verbose) {
      message("# Segmenting Airways: Identifying Airway Location")
    }
    fimg <- air_mask > 2
    mydim = dim(fimg)
    fimg<-as.array(fimg,dim=mydim)
    fimg[round(mydim[1]/4):(round(mydim[1]/4)*3),round(mydim[2]/4):(round(mydim[2]/4)*3),round(mydim[3]/2):mydim[3]] <- 1
    fimg<-makeImage(mydim,fimg)
    air_mask <- maskImage(air_mask,fimg)

    if (verbose) {
      message("# Segmenting Airways: Smoothing")
    }
    air_mask = iMath(air_mask, "MC", 1)
    air_mask = iMath(air_mask, "ME", 1)
    air_mask = iMath(air_mask, "MD", 2)


    # Segmenting Lung
    if (verbose) {
      message("# Segmenting Lung: Removing Airways from Lung")
    }
    lung_mask = maskImage(lung_air_mask, 1-air_mask)
    img_masked = maskImage(img, lung_mask)

    if (verbose) {
      message("# Segmenting Lung: Finding Left and Right Lungs")
    }
    left_right_mask = labelClusters(lung_mask, minClusterSize = 50000)
    n_clus = length(unique(left_right_mask))

    if (n_clus == 1) {
      message("# Error: Can't find lungs, returning current mask")
      return(lung_mask)
    }

    # Left and right lung are still connected
    i = 0
    j = 0
    lung_mask2 = antsImageClone(lung_mask)
    while(n_clus == 2){
      if(i == 0) {
        lung_mask2 = iMath(lung_mask2, "ME", 1)
        j = j + 1
        left_right_mask = labelClusters(lung_mask2, minClusterSize = 50000)
        n_clus = length(unique(left_right_mask))
      }
      i = i + 1
      new_lthresh = lthresh - 50*i
      if(new_lthresh < 0){
        message("# Error: Can't distinguish left/right lungs, returning left/right combined mask")
        return(lung_mask)
      }
      if(new_lthresh <= 200) {
        lung_mask2 = iMath(lung_mask2, "ME", 1)
        j = j + 1
      }
      img_mask = img_masked < new_lthresh
      lung_mask2[img_mask == 0] = 0
      left_right_mask = labelClusters(lung_mask2, minClusterSize = 50000)
      num = unique(left_right_mask)
      n_clus = length(num)
    }

    # Correct number of clusters
    if (n_clus == 3) {

      # Find which value is the background (Note: this was added because sometimes the background is non-zero)
      left_right_mask = left_right_mask + 1
      fimg <- as.array(left_right_mask)
      last <- dim(fimg)[1]
      left <- unique(as.vector(fimg[1,,]))
      right <- unique(as.vector(fimg[last,,]))
      both <- c(left,right)
      # Logic: only the background value will be on both the left and right sides
      value <- both[duplicated(both)]

      # Make background value 0, and the other two values 1 and 2
      left_right_mask[left_right_mask == value] <- 0
      num = unique(left_right_mask)
      if(1 %in% num & 3 %in% num){
        left_right_mask[left_right_mask==3] <- 2
      }
      if(2 %in% num & 3 %in% num){
        left_right_mask[left_right_mask==2] <- 1
        left_right_mask[left_right_mask==3] <- 2
      }


      # Find coordinates of left and right clusters
      coord = sapply(1:2, function(i){
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
        message("# Segmenting Lung: Finishing Touches")
      }
      left_mask = iMath(left_mask, "MC", 8+2*j)
      left_mask = iMath(img = left_mask, operation = "GetLargestComponent")
      right_mask = iMath(right_mask, "MC", 8+2*j)
      right_mask = iMath(img = right_mask, operation = "GetLargestComponent")
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

























