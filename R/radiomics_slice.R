#' Radiomic Calculation on CT Slices
#'
#' Calculate radiomic features on the each 2D slice of the whole 3D lung, right and left lungs separately
#'
#' @param img CT scan in ANTs image file format
#' @param mask Mask of CT scan in ANTs image file format
#' @param sides Choose to calculate radiomic features on the right and/or left lungs. Note: Right lung = 1, left lung = 2, non-lung = 0
#' @param plane One of: axial, coronal, sagittal
#' @param featuresFirst First level radiomic features to calculate
#' @param featuresSpatial Spatial radiomic features to calculate
#' @param tidy Logical. If true, outputs a tidy dataframe with results. If false, outputs nested loop.
#' @param reduce Logical. If true, reduces the dimensions of the scan based on extent of mask using reduce_scan.
#'
#' @return Radiomic values from every slice in both the right and left lungs
#' @export
#'
#' @examples
radiomics_slice <- function(img,
                           mask,
                           sides = c("right", "left"),
                           plane = 'axial',
                           featuresFirst = c('mean', 'sd', 'skew', 'kurtosis', 'min', 'q1', 'median', 'q3', 'max','energy', 'rms', 'uniformity', 'entropy'),
                           featuresSpatial = c('mi', 'gc', 'fd'),
                           tidy = TRUE,
                           reduce = TRUE){


  featuresMask <- lapply(sides, function(side){

    if(side == "right"){mv = 1}
    if(side == "left"){mv = 2}

    mask2 <- mask == mv

    # Reduce scan (optional)
    if(reduce == TRUE){
      red <- reduce_scan(img, mask2)
      img2 <- red$img
      mask2 <- red$mask
      rm(red)
      gc()
    }else{img2 <- img}

    # Put image in array format and remove non-mask values
    img2 <- as.array(img2)
    mask2 <- as.array(mask2)
    img2[mask2 != 1] <- NA


    # Which plane?
    if(plane == 'axial'){p = 3}
    if(plane == 'coronal'){p = 2}
    if(plane == 'sagittal'){p = 1}


    # Calculate features
    ndim <- dim(img2)[p]
    features <- apply(img2, p, function(x){
      npixels <- length(x[!is.na(x)])
      if(length(featuresFirst)>0){
        features1 <- radiomics_first(x, featuresFirst)
      }else(features1 <- NULL)
      if(length(featuresSpatial)>0){
        features2 <- radiomics_spatial(x, featuresSpatial)
      }else(features2 <- NULL)
      features <- c(features1, features2)
      features <- features[c(featuresFirst, featuresSpatial)]
      features <- c(npixels = npixels, features)
      return(features)
    })
    names(features) <- paste0('slic_num_', 1:ndim)

    return(features)
  })
  names(featuresMask) <- sides


  if(tidy == TRUE){
    # Make a nice little data frame to output
    test2 = NULL
    for(i in 1:length(featuresMask)){
      test <- do.call('rbind', featuresMask[[i]])
      test <- cbind.data.frame(lung = names(featuresMask)[i],
                               slice_number = names(featuresMask[[i]]),
                               test)
      test$slice_number <- gsub("slic_num_", "", test$slice_number)
      test <- as.data.frame(sapply(test, as.numeric))
      nslic <- dim(test)[1]
      test$slice_percent <- test$slice_number/nslic * 100
      test2 <- rbind(test2, test)
    }
    featuresMask <- test2
    rownames(featuresMask) <- c()
  }

  return(featuresMask)
}