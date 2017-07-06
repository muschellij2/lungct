#' @title Segment Human from CT scans
#' @description Segment Human from a non-contrast CT scan
#'
#' @param img Filename or \code{antsImage}
#' @param adder amount to be added to the image to make non-zero
#' @param lthresh lower threshold for the image
#' @param verbose Print diagnostic messages
#'
#' @return List of smoothed image, body, adder
#' @export
segment_human = function(
  img,
  adder = 1025,
  lthresh = -300,
  verbose = TRUE
) {
  reg_img = check_ants(img)
  ##############################
  # 1024 should be lower limit
  ##############################
  if (verbose) {
    message("Making Positive Values")
  }
  # reg_img = as.array(reg_img)

  reg_img = reg_img + adder
  if (verbose) {
    message("Setting voxel ranges")
  }
  reg_img[reg_img < 0] = 0
  reg_img[reg_img > 3071 + adder] = 3071 + adder
  # reg_img = as.antsImage(reg_img, reference = img)

  if (verbose) {
    message("# Getting Humans: Smoothing Image")
  }
  ss = iMath(reg_img, "PeronaMalik", 10, 5)

  if (verbose) {
    message("# Getting Humans: Largest Component")
  }
  body = ss > (0 + adder)
  body = iMath(img = body, operation = "GetLargestComponent")
  inds = getEmptyImageDimensions(body)
  if (verbose) {
    message("# Getting Humans: Dropping Zero Dimensions")
  }
  ss = maskEmptyImageDimensions(
    img = ss,
    inds = inds,
    mask.value = 0)

  if (verbose) {
    message("# Getting Humans: Making Coarse Body")
  }
  body = ss > (lthresh + adder)
  body = iMath(img = body, operation = "GetLargestComponent")

  L = list(
    smoothed = ss,
    body = body,
    adder = adder)
  return(L)
}