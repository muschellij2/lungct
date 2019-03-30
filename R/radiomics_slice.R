#' Calculate radiomic features on the each 2D slice of the whole 3D lung, left and right lungs separately
#'
#' @param img CT scan in ANTs image file format
#' @param mask Mask of CT scan in ANTs image file format
#' @param background_value Value of background
#' @param plane One of: axial, coronal, sagittal
#' @param featuresFirst First level radiomic features to calculate
#' @param featuresSpatial Spatial radiomic features to calculate
#' @param tidy Logical. If true, outputs a tidy dataframe with results. If false, outputs nested loop.
#'
#' @return Radiomic values from every slice in both the left and right lungs
#' @importFrom data.table rbindlist
#' @export
#'
#' @examples
radiomics_slice <- function(img,
                           mask,
                           background_value = NA,
                           plane = 'axial',
                           featuresFirst = c('mean', 'sd', 'skew', 'kurtosis', 'min', 'q1', 'median', 'q3', 'max','energy', 'rms', 'uniformity', 'entropy'),
                           featuresSpatial = c('mi', 'gc', 'fd'),
                           tidy = TRUE){

  # Find unique mask values
  mask_values <- unique(mask)
  mask_values <- mask_values[-which(mask_values == background_value)]


  featuresMask <- lapply(mask_values, function(mv){

    # Put image in array format and remove non-mask values
    img2 <- as.array(img)
    mask2 <- as.array(mask)
    img2[mask2 != mv] <- NA


    # Which plane?
    if(plane == 'axial'){p = 3}
    if(plane == 'coronal'){p = 2}
    if(plane == 'sagittal'){p = 1}


    # Calculate features
    ndim <- dim(img2)[p]
    features <- apply(img2, p, function(x){
      npixels <- length(x[!is.na(x)])
      features1 <- radiomics_first(x, featuresFirst)
      features2 <- radiomics_spatial(x, featuresSpatial)
      features <- c(npixels = npixels, features1, features2)
      return(features)
    })
    names(features) <- paste0('slic_num_', 1:ndim)

    return(features)
  })
  names(featuresMask) <- paste0('mask',mask_values)

  if(tidy == TRUE){
    # Make a nice little data frame to output
    test2 = NULL
    for(i in 1:length(featuresMask)){
      test <- rbindlist(featuresMask[[i]])
      test <- cbind.data.frame(mask = names(featuresMask)[i],
                               slic_num = names(featuresMask[[i]]),
                               test)
      test2 <- rbind(test2, test)
    }
    featuresMask <- test2
  }

  return(featuresMask)
}