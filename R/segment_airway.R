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
  # air_mask = labelClusters(air_mask, minClusterSize = 5000)
  # n_clus = length(unique(air_mask))
  #
  # if(n_clus > 2){
  #   coord = sapply(1:n_clus, function(i){
  #     iimg = air_mask == i
  #     loc = which(as.array(iimg) == 1, arr.ind = T)
  #     x = median(loc[,1])
  #     return(x)
  #   })
  #   center = dim(air_mask)[1]/2
  #   loss = abs(center - coord)
  #   clus = which(loss == min(loss, na.rm = T))
  #   air_mask = air_mask == clus
  # }
  air_mask = iMath(air_mask, "MD", 2)

  return(air_mask)
}
