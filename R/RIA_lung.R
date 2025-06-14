#' Calculate first-order, GLCM, and/or GLRLM radiomic features on the whole 3D lung
#'
#' This is a wrapper for \code{RIA} R package. It calculates first-order, GLCM, and/or GLRLM on the whole 3D lung, left and right lungs separately
#'
#' @param img CT scan in ANTs image file format
#' @param mask Mask of CT scan in ANTs image file format
#' @param sides Choose to calculate radiomic features on the right and/or left lungs. Note: Right lung = 1, left lung = 2, non-lung = 0
#' @param features Type of radiomic features to calculate. Options: first-order, GLCM, and/or GLRLM
#' @param bins_in Number of bins to discretize image
#' @param equal_prob logical, indicating to cut data into bins with equal relative frequencies.
#' If FALSE, then equal interval bins will be used.
#' @param distance integer, distance between the voxels being compared.
#' @param statistic string, defining the statistic to be calculated on the array of GLCM statistics.
#' By default, statistic is set to \emph{"mean"}, however any function may be provided. The proper
#' syntax is: function(X, attributes). The supplied string must contain a "X", which will be replaced
#' with the array of the GLCM statistics value. Further attributes of the function may also be given.
#' For example, if you wish to calculate the median of all GLCMs calculated in different directions,
#' then it must be supplied as: \emph{median(X, na.rm = TRUE)}.
#' @param verbose_in logical, indicating whether to print detailed information.
#' Most prints can also be suppressed using the \code{\link{suppressMessages}} function.
#'
#' @return list containing the statistical information
#' @export
RIA_lung <- function(img,
                     mask,
                     sides = c("right", "left"),
                     features = c('fo', 'glcm', 'glrlm'),
                     bins_in = 8,
                     equal_prob = FALSE,
                     distance = 1,
                     statistic = "mean(X, na.rm = TRUE)",
                     verbose_in = TRUE){

  if (!requireNamespace("RIA", quietly = TRUE)) {
    stop("RIA package required for RIA_lung")
  }
  # Loop through mask values
  featuresMask <- lapply(sides, function(side){

    if(side == "right"){mv = 1}
    if(side == "left"){mv = 2}

    # Put image in array format and remove non-mask values
    data <- as.array(img)
    mask2 <- as.array(mask)
    data[mask2 != mv] <- NA

    no_na_sum = function(x) {
      sum(x, na.rm = TRUE)
    }
    # Crop image to speed up computation
    test <- apply(data, 1, no_na_sum)
    data <- data[which(test != 0),,]
    test <- apply(data, 2, no_na_sum)
    data <- data[,which(test != 0),]
    test <- apply(data, 3, no_na_sum)
    data <- data[,,which(test != 0)]

    ###create RIA_image structure
    RIA_image <- list(data = NULL, header = list(), log = list())
    if(length(dim(data)) == 3 | length(dim(data)) == 2) {class(RIA_image) <- append(class(RIA_image), "RIA_image")
    } else {
      stop(paste0("ANTsImage LOADED IS ", length(dim(data)),
                  " DIMENSIONAL. ONLY 2D AND 3D DATA ARE SUPPORTED!"))
    }
    RIA_image$data$orig  <- data
    RIA_image$data$modif <- NULL
    class(RIA_image$header) <- append(class(RIA_image$header), "RIA_header")
    class(RIA_image$data) <- append(class(RIA_image$data), "RIA_data")
    class(RIA_image$log) <- append(class(RIA_image$log), "RIA_log")
    RIA_image$log$events  <- "Created"
    RIA_image$log$orig_dim  <- dim(data)


    # Calculate first order radiomic features
    if('fo' %in% features){
      RIA_image <- RIA:::first_order(RIA_image, use_type = "single", use_orig = TRUE, verbose_in = verbose_in)
    }


    # Discretize image
    if('glcm' %in% features | 'glrlm' %in% features){
      RIA_image <- RIA:::discretize(RIA_image, bins_in=bins_in, equal_prob = equal_prob, verbose_in = verbose_in)

      # Calculate GLCM radiomic features
      if('glcm' %in% features){
        for (i in 1: length(distance)) {
          RIA_image <- RIA:::glcm_all(RIA_image, use_type = "discretized", distance = distance[i], verbose_in = verbose_in)
        }
        RIA_image <- RIA:::glcm_stat(RIA_image, use_type = "glcm", verbose_in = verbose_in)
        RIA_image <- RIA:::glcm_stat_all(RIA_image, statistic = statistic, verbose_in = verbose_in)
      }


      # Calculate GLRLM radiomic features
      if ('glrlm' %in% features){
        RIA_image <- RIA:::glrlm_all(RIA_image, use_type = "discretized", verbose_in = verbose_in)
        RIA_image <- RIA:::glrlm_stat(RIA_image, use_type = "glrlm", verbose_in = verbose_in)
        RIA_image <- RIA:::glrlm_stat_all(RIA_image, statistic = statistic, verbose_in = verbose_in)
      } else{
        RIA_image$stat_glrlm_mean <- NULL
      }
    }


    features <- list(first_order = RIA_image$stat_fo$orig,
                     glcm = RIA_image$stat_glcm_mean,
                     glrlm = RIA_image$stat_glrlm_mean)


    return(features)
  })
  names(featuresMask) <- sides

  return(featuresMask)
}