#' @title Segment Lungs from CT scans
#' @description Segment Lungs from a non-contrast CT scan
#'
#' @param img Filename or \code{antsImage}
#' @param verbose Print diagnostic messages
#'
#' @return List of image, lung, lung mask
#' @export
#' @importFrom extrantsr check_ants filler
#' @importFrom neurobase zero_pad maskEmptyImageDimensions
#' @importFrom oro.nifti voxdim
#' @importFrom stats median
#' @importFrom ANTsRCore resampleImage iMath smoothImage as.antsImage
#' @importFrom ANTsRCore antsImageClone
segment_lung = function(img, verbose = TRUE) {

  reg_img = check_ants(img)
  img = antsImageClone(reg_img)
  vres = voxdim(reg_img)
  reg_img = resampleImage(reg_img, c(1,1,1))


  ##############################
  # 1024 should be lower limit
  ##############################
  adder = 1025
  reg_img = reg_img + adder
  reg_img[reg_img < 0] = 0
  reg_img[reg_img > 3071 + adder] = 3071 + adder

  if (verbose) {
    message("# Getting Humans")
  }
  ss = iMath(reg_img, "PeronaMalik", 10, 5)
  body = ss > (0 + adder)
  body = iMath(img = body, operation = "GetLargestComponent")
  inds = getEmptyImageDimensions(body)
  rm(list = "body"); gc(); gc()
  ss = maskEmptyImageDimensions(
    img = ss,
    inds = inds,
    mask.value = 0)

  if (verbose) {
    message("# smoothing image")
  }
  ss = smoothImage(inimg = ss, sigma = 10,
                   max_kernel_width = 200)
  # smoothed image is greater than some thresh

  med = median(as.numeric(ss))
  body = ss > med

  # body = ss > (-100 + adder)
  if (verbose) {
    message("# Getting Humans")
  }

  # zero padding in case for any connectedness
  zp = as.array(body)
  kdim = c(1,1,1)
  zp = zero_pad(zp, kdim = kdim)
  zp = as.antsImage(zp, reference = body)
  # getting connected component
  cc = iMath(img = zp, operation = "GetLargestComponent")
  rm(list = "zp"); gc(); gc()
  cc = as.array(cc)
  cc = zero_pad(cc, kdim = kdim, invert = TRUE)
  cc = as.antsImage(cc, reference = body)

  if (verbose) {
    message("# Filling Holes")
  }
  # cc = iMath(img = cc, operation = "FillHoles")
  cc = filler(cc, fill_size = 60)
  # cc = filler(cc, fill_size = 40)
  cc = filler(cc, fill_size = 5, dilate = FALSE)

  # Dropping non-human stuff
  inds = getEmptyImageDimensions(cc)

  if (verbose) {
    message("# Making New image")
  }

  # Don't want to drop the indices,
  # just blank them out
  newimg = maskEmptyImageDimensions(img = reg_img, inds = inds)

  lung = newimg < (-300 + adder) & newimg > 0 & cc == 1
  lung = iMath(img = lung, operation = "GetLargestComponent")
  # lung = iMath(img = lung, operation = "FillHoles")
  lung = filler(lung, fill_size = 2)
  lung_mask = resampleImage(lung,
                            resampleParams = dim(img), useVoxels = TRUE,
                            interpType = 1)
  # lung_mask = resampleImageToTarget(
  #   lung, target = img,
  #   interpType = "nearestNeighbor",
  #   verbose = verbose)
  lung = maskImage(img, lung_mask)

  # reg_img = reg_img - adder
  L = list(img = img,
           lung_mask = lung_mask,
           lung = lung
  )
  return(L)
}
