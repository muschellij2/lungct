#' Calculate first order radiomic features on a 2D or 3D array
#'
#' @param data Any 2D or 3D image (as matrix or array) to calculate first-order features
#' @param features = first level radiomic features to calculate
#'
#' @return Values from selected features
#' @importFrom stats quantile sd
#' @export
radiomics_first <- function(data,
                       features = c('mean', 'sd', 'skew', 'kurtosis', 'min', 'q1', 'median', 'q3', 'max','energy', 'rms', 'uniformity', 'entropy')){

  # Clean up data
  data <- as.vector(data)
  data <- data[!is.na(data)]

  # Moments
  if('mean' %in% features){
    mean_value <- mean(data)
  }else(mean_value = NULL)
  if('sd' %in% features){
    sd_value <- sd(data)
  }else(sd_value = NULL)
  if('skew' %in% features){
    skew_value <- skew(data)
  }else(skew_value = NULL)
  if('kurtosis' %in% features){
    kurtosis_value <- kurtosis(data)
  }else(kurtosis_value = NULL)


  # Quartiles
  if('min' %in% features){
    min_value <- min(data)
  }else(min_value = NULL)
  if('q1' %in% features){
    q1_value <- unname(quantile(data, probs = 0.25))
  }else(q1_value = NULL)
  if('median' %in% features){
    median_value <- median(data)
  }else(median_value = NULL)
  if('q3' %in% features){
    q3_value <- unname(quantile(data, probs = 0.75))
  }else(q3_value = NULL)
  if('max' %in% features){
    max_value <- max(data)
  }else(max_value = NULL)


  # Others
  if('energy' %in% features){
    energy_value <- energy(data)
  }else(energy_value = NULL)
  if('rms' %in% features){
    rms_value <- rms(data)
  }else(rms_value = NULL)
  if('uniformity' %in% features){
    uniformity_value <- uniformity(data)
  }else(uniformity_value = NULL)
  if('entropy' %in% features){
    entropy_value <- entropy(data)
  }else(entropy_value = NULL)


  featuresList <- list(
    mean = mean_value,
    sd = sd_value,
    skew = skew_value,
    kurtosis = kurtosis_value,
    min = min_value,
    q1 = q1_value,
    median = median_value,
    q3 = q3_value,
    max = max_value,
    energy = energy_value,
    rms = rms_value,
    uniformity = uniformity_value,
    entropy = entropy_value
  )
  if(length(features)==1){
    featuresList = unlist(featuresList[features])
  }

  return(featuresList)

}



skew <- function(data) {
  avg <- mean(data)
  SD  <- stats::sd(data)
  output <- mean(((data-avg)^3))/(SD)^3
  return(output)
}

kurtosis <- function(data) {
  avg <- mean(data)
  SD  <- stats::sd(data)
  output <- mean(((data-avg)^4))/(SD)^4
  return(output-3)
}

energy <- function(data) {
  output <- sum(data^2)
  return(output)
}

rms <- function(data) {
  output <- sqrt(sum(data^2)/length(data))
  return(output)
}

uniformity <- function(data) {
  output <- sum((table(data)/length(data))^2)
  return(output)
}

entropy <- function (data, base = 2)
{
  p <- table(data)/length(data)
  l <- logb(p, base)
  H <- sum(p * l)*-1
  return(H)
}




















