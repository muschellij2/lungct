#' Calculate spatial radiomic features on a 2D or 3D array
#'
#' @param data Any 2D or 3D image (as matrix or array) to calculate spatial features
#' @param features = spatial radiomic features to calculate
#' @return Values from selected features
#' @importFrom abind abind
#' @export
radiomics_spatial <- function(data,
                              features = c('mi', 'gc')){


  # Figure data dimension
  dimData <- length(dim(data))

  #2D
  if(dimData == 2){
    if('mi' %in% features){
      mi_value <- moran2D(data)
    }else(mi_value = NULL)
    if('gc' %in% features){
      gc_value <- geary2D(data)
    }else(gc_value = NULL)

  }

  # 3D
  if(dimData == 3){

    if('mi' %in% features){
      mi_value <- moran3D(data)
    }else(mi_value = NULL)
    if('gc' %in% features){
      gc_value <- geary3D(data)
    }else(gc_value = NULL)

  }

  features <- list(
    mi = mi_value,
    gc = gc_value
  )

  return(features)

}


# Calculate Moran's I in 3D
moran3D <- function(mat){

  # Set up
  n <- length(mat[is.na(mat)==F])
  xbar <- mean(mat, na.rm=T)
  diff <- mat - xbar
  diff2 <- diff**2
  sum_diff2 <- sum(diff2, na.rm = T)

  # Dimension of matrix
  nx <- dim(diff)[1]
  ny <- dim(diff)[2]
  nz <- dim(diff)[3]

  # Rook shift
  left <- diff * abind(diff[-1,,], matrix(NA, nrow = ny, ncol = nz), along = 1)
  right <- diff * abind(matrix(NA, nrow = ny, ncol = nz), diff[-nx,,], along = 1)
  front <- diff * abind(diff[,-1,], matrix(NA, nrow = nx, ncol = nz), along = 2)
  back <- diff * abind(matrix(NA, nrow = nx, ncol = nz), diff[,-ny,], along = 2)
  up <- diff * abind(diff[,,-1], matrix(NA, nrow = nx, ncol = ny), along = 3)
  down <- diff * abind(matrix(NA, nrow = nx, ncol = ny), diff[,,-nz], along = 3)

  # Calcultate weights
  weights <- sum(length(left[is.na(left) == F])) +
    sum(length(right[is.na(right) == F])) +
    sum(length(front[is.na(front) == F])) +
    sum(length(back[is.na(back) == F])) +
    sum(length(up[is.na(up) == F])) +
    sum(length(down[is.na(down) == F]))

  # Change all "NAs" to zero, so we can efficiently sum matrices
  left[is.na(left)] <- 0
  right[is.na(right)] <- 0
  front[is.na(front)] <- 0
  back[is.na(back)] <- 0
  up[is.na(up)] <- 0
  down[is.na(down)] <- 0

  # Put info together to calculate MI
  sumv <- left + right + front + back + up + down
  local_mi <- sumv / sum_diff2 * n / weights
  global_mi <- sum(local_mi)

  return(global_mi)
}




# Calculate Moran's I in 2D
moran2D<-function(mat){

  # Set up
  n <- length(mat[is.na(mat)==F])
  xbar <- mean(mat, na.rm=T)
  diff <- mat - xbar
  diff2 <- diff**2
  sum_diff2 <- sum(diff2, na.rm = T)

  # Dimension of matrix
  nrow <- dim(diff)[1]
  ncol <- dim(diff)[2]

  # Rook shift
  up <- diff * rbind(diff[-1,],rep(NA,ncol))
  down <- diff * rbind(rep(NA,ncol),diff[-nrow,])
  left <- diff * cbind(diff[,-1],rep(NA,nrow))
  right <- diff * cbind(rep(NA,nrow),diff[,-ncol])

  # Calcultate weights
  weights <- sum(length(left[is.na(left) == F])) +
    sum(length(right[is.na(right) == F])) +
    sum(length(up[is.na(up) == F])) +
    sum(length(down[is.na(down) == F]))

  # Change all "NAs" to zero, so we can efficiently sum matrices
  left[is.na(left)] <- 0
  right[is.na(right)] <- 0
  up[is.na(up)] <- 0
  down[is.na(down)] <- 0

  # Put info together to calculate MI
  sumv <- left + right + up + down
  local_mi <- sumv / sum_diff2 * n / weights
  global_mi <- sum(local_mi)

  return(global_mi)
}




# Calculate Geary's C in 3D
geary3D<-function(mat){

  # Set up
  n <- length(mat[is.na(mat)==F])
  xbar <- mean(mat, na.rm=T)
  diff <- mat - xbar
  diff2 <- diff**2
  sum_diff2 <- sum(diff2, na.rm = T)

  # Dimension of matrix
  nx <- dim(diff)[1]
  ny <- dim(diff)[2]
  nz <- dim(diff)[3]

  # Rook shift
  left <- (mat - abind(mat[-1,,], matrix(NA, nrow = ny, ncol = nz), along = 1) )^2
  right <- (mat - abind(matrix(NA, nrow = ny, ncol = nz), mat[-nx,,], along = 1) )^2
  front <- (mat - abind(mat[,-1,], matrix(NA, nrow = nx, ncol = nz), along = 2) )^2
  back <- (mat - abind(matrix(NA, nrow = nx, ncol = nz), mat[,-ny,], along = 2) )^2
  up <- (mat - abind(mat[,,-1], matrix(NA, nrow = nx, ncol = ny), along = 3) )^2
  down <- (mat - abind(matrix(NA, nrow = nx, ncol = ny), mat[,,-nz], along = 3) )^2

  # Calcultate weights
  weights <- sum(length(left[is.na(left) == F])) +
    sum(length(right[is.na(right) == F])) +
    sum(length(front[is.na(front) == F])) +
    sum(length(back[is.na(back) == F])) +
    sum(length(up[is.na(up) == F])) +
    sum(length(down[is.na(down) == F]))

  # Change all "NAs" to zero, so we can efficiently sum matrices
  left[is.na(left)] <- 0
  right[is.na(right)] <- 0
  front[is.na(front)] <- 0
  back[is.na(back)] <- 0
  up[is.na(up)] <- 0
  down[is.na(down)] <- 0

  # Put info together to calculate GC
  sumv <- left + right + front + back + up + down
  local_gc <- sumv / sum_diff2 * (n - 1) / (2 * weights)
  global_gc <- sum(local_gc)

  return(global_gc)
}


# Calculate Geary's C in 2D
geary2D<-function(mat){

  # Set up
  n <- length(mat[is.na(mat)==F])
  xbar <- mean(mat, na.rm=T)
  diff <- mat - xbar
  diff2 <- diff**2
  sum_diff2 <- sum(diff2, na.rm = T)

  # Dimension of matrix
  nrow <- dim(diff)[1]
  ncol <- dim(diff)[2]

  # Rook shift
  up <- (mat - rbind(mat[-1,],rep(NA,ncol)) )^2
  down <- (mat - rbind(rep(NA,ncol),mat[-nrow,]) )^2
  left <- (mat - cbind(mat[,-1],rep(NA,nrow)) )^2
  right <- (mat - cbind(rep(NA,nrow),mat[,-ncol]) )^2

  # Calcultate weights
  weights <- sum(length(left[is.na(left) == F])) +
    sum(length(right[is.na(right) == F])) +
    sum(length(up[is.na(up) == F])) +
    sum(length(down[is.na(down) == F]))

  # Change all "NAs" to zero, so we can efficiently sum matrices
  left[is.na(left)] <- 0
  right[is.na(right)] <- 0
  up[is.na(up)] <- 0
  down[is.na(down)] <- 0

  # Put info together to calculate GC
  sumv <- left + right + up + down
  local_gc <- sumv / sum_diff2 * (n - 1) / (2 * weights)
  global_gc <- sum(local_gc)

  return(global_gc)
}



















