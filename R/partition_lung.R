#' Partition the lung
#'
#' @param img Lung scan to be partitioned. Can be in 3D matrix or ANTs image file format.
#' @param kernel_size Size of the kernel, in voxel units of width, depth, and height. Must be c(3,3,3) or greater. Default: c(30,30,30)
#' @param kernel_stride Stride (or spacing) between kernels, in voxel units, for width, depth, and height. If kernel_stride = kernel_size, the partitions are non-overlapping. If stride = c(1,1,1), then each voxel is returned.
#' @param centroid Logical. If true, output includes the centroids of each partition.
#'
#' @return Matrix of x, y, and z coordinates for each partition
#' @export
partition_lung = function(img,
                          kernel_size = c(30,30,30),
                          kernel_stride = c(30,30,30),
                          centroid = TRUE){

  # Put in array format
  img <- as.array(img)

  # Dimensions of image
  dim_x <- dim(img)[1]
  dim_y <- dim(img)[2]
  dim_z <- dim(img)[3]

  # Kernel sequences -- start points
  x1 = seq(1, dim_x, by = kernel_stride[1])
  y1 = seq(1, dim_y, by = kernel_stride[2])
  z1 = seq(1, dim_z, by = kernel_stride[3])

  # Kernel endpoints
  xend <- x1 + kernel_size[1] - 1
  yend <- y1 + kernel_size[2] - 1
  zend <- z1 + kernel_size[3] - 1

  # Kernel endpoints can't be larger than dimension of image
  xend[xend > dim_x] = dim_x
  yend[yend > dim_y] = dim_y
  zend[zend > dim_z] = dim_z

  # Put into dataframe
  coords1 <- expand.grid(x1, y1, z1)
  coords_end <- expand.grid(xend, yend, zend)
  coords1 <- cbind(coords1, coords_end)
  colnames(coords1) <- c('x1','y1','z1','xend','yend','zend')
  coords1$partition <- 1:dim(coords1)[1]

  # Find the centroid for each partition
  if(centroid == TRUE){
    coords1$x <- apply(cbind(coords1$x1, coords1$xend), 1, median)
    coords1$y <- apply(cbind(coords1$y1, coords1$yend), 1, median)
    coords1$z <- apply(cbind(coords1$z1, coords1$zend), 1, median)
  }

  # Return coordinates
  return(coords1)

}
