#' Calculate radiomic features on the partitioned 3D lung
#'
#' @param img CT scan in ANTs image file format
#' @param mask Mask of CT scan in ANTs image file format
#' @param mask_values Values of mask to use for radiomic feature calculation
#' @param featuresFirst First level radiomic features to calculate
#' @param featuresSpatial Spatial radiomic features to calculate
#' @param partition Matrix of x, y, and z coordinates for each partition from partition_lung. If null, partition_lung is called.
#' @param kernel_size (If partition is null) Size of the kernel, in voxel units of width, depth, and height. Must be c(3,3,3) or greater. Default: c(30,30,30)
#' @param kernel_stride (If partition is null) Stride (or spacing) between kernels, in voxel units, for width, depth, and height. If kernel_stride = kernel_size, the partitions are non-overlapping. If stride = c(1,1,1), then each voxel is returned.
#' @param threshold Number of non-missing voxels needed to calculate radiomic features in each partition.
#' @param tidy Logical. If true, outputs a tidy dataframe with results. If false, outputs nested loop.
#'
#' @return Values from selected features for both left and right lungs
#' @importFrom ANTsR maskImage
#' @export
radiomics_partition <- function(img,
                                mask,
                                mask_values = c(1, 2),
                                featuresFirst = c('mean', 'sd', 'skew', 'kurtosis', 'min', 'q1', 'median', 'q3', 'max','energy', 'rms', 'uniformity', 'entropy'),
                                featuresSpatial = c('mi', 'gc', 'fd'),
                                partition = NULL,
                                kernel_size = c(30, 30, 30),
                                kernel_stride = c(30, 30, 30),
                                threshold = 1000,
                                tidy = TRUE) {

  # Get partition, if necessary
  if(is.null(partition)){
    partition = partition_lung(img,
                               kernel_size = kernel_size,
                               kernel_stride = kernel_stride,
                               centroid = TRUE)
  }


  # Calculate radiomic features on partitions within each mask value
  featuresMask <- lapply(mask_values, function(mv){


    # Put image in array format and remove non-mask values
    img2 <- as.array(img)
    mask2 <- as.array(mask)
    mask2 <- mask2 == mv
    img2[mask2 != 1] <- NA

    # Calculate n each partition
    features <- lapply(1:dim(partition)[1], function(i){

      # Grab partition
      x <- img2[partition$x1[i]:partition$xend[i],
                partition$y1[i]:partition$yend[i],
                partition$z1[i]:partition$zend[i]]

      # Find dimension of partition and number of non-null pixels
      dim_p <- dim(x)
      if(is.null(dim_p)){dim_p <- c(0,0,0)}
      npixels <- length(x[!is.na(x)])


      # Only calculate radiomic features if partition fits criteria
      if(dim_p[1] > 2 & dim_p[2] > 2 & dim_p[3] > 2 & npixels >= threshold){

        # Calculate features
        if(length(featuresFirst)>0){
          features1 <- radiomics_first(x, featuresFirst)
        }else(features1 <- NULL)
        if(length(featuresSpatial)>0){
          features2 <- radiomics_spatial(x, featuresSpatial)
        }else(features2 <- NULL)


        # Put features together and only keep specified features
        features <- c(features1, features2)
        features <- features[c(featuresFirst, featuresSpatial)]
        features <- c(npixels = npixels, features)

      }else(features <- NULL)

      return(features)
    })

    # Name partitions and remove NULL partitions
    names(features) <- paste0('partition', 1:dim(partition)[1])
    # features[sapply(features, is.null)] <- NULL
    return(features)
  })
  names(featuresMask) <- paste0('mask',mask_values)


  if(tidy == TRUE){
    # Make a nice little data frame to output
    test2 = NULL
    for(i in 1:length(featuresMask)){

      # Reduce list to data frame
      test <- do.call('rbind', featuresMask[[i]])
      test <- cbind.data.frame(mask_value = names(featuresMask)[i],
                               partition = names(featuresMask[[i]]),
                               test)

      # Reformatting
      test$partition <- gsub("partition", "", test$partition)
      test$mask_value <- gsub("mask", "", test$mask_value)
      test <- as.data.frame(sapply(test, as.numeric))

      # Add in partition centroids
      partition2 <- partition[partition$partition %in% test$partition, 8:10]
      test <- cbind(partition2, test)
      test2 <- rbind(test2, test)
    }
    featuresMask <- test2
    rm(test,test2,partition2)
    gc()
    rownames(featuresMask) <- c()
  }

  return(featuresMask)
}



















