#' Calculate individual radiomic features on the whole 3D lung, left and right lungs separately
#'
#' @param img CT scan in ANTs image file format
#' @param mask Mask of CT scan in ANTs image file format
#' @param background_value = value of background
#' @param featuresFirst = first level radiomic features to calculate
#' @param featuresSpatial = spatial radiomic features to calculate
#'
#' @return Values from selected features for both left and right lungs
#' @importFrom ANTsR maskImage
#' @export
radiomics_lung <- function(img,
                          mask,
                          background_value = NA,
                          featuresFirst = c('mean', 'sd', 'skew', 'kurtosis', 'min', 'q1', 'median', 'q3', 'max','energy', 'rms', 'uniformity', 'entropy'),
                          featuresSpatial = c('mi', 'gc')){

  # Find unique mask values
  mask_values <- unique(mask)
  mask_values <- mask_values[-which(mask_values == background_value)]


  featuresMask <- lapply(mask_values, function(mv){

    # Put image in array format and remove non-mask values
    img2 <- as.array(img)
    mask2 <- as.array(mask)
    img2[mask2 != mv] <- NA

    # Calculate features
    features1 <- radiomics_first(img2, featuresFirst)
    features2 <- radiomics_spatial(img2, featuresSpatial)
    return(c(features1, features2))
  })
  names(featuresMask) <- paste0('mask',mask_values)

  return(featuresMask)
}



















