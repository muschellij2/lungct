#' @title Segment Lungs from CT scans
#' @description Segment Lungs from a non-contrast CT scan
#'
#' @param img Filename or \code{antsImage}
#' @param mask Should the image be masked
#' @param lthresh lower threshold for the image
#' @param verbose Print diagnostic messages
#'
#' @return List of image, lung, lung mask
#' @export
#' @importFrom extrantsr check_ants filler
#' @importFrom neurobase zero_pad maskEmptyImageDimensions mask_img
#' @importFrom oro.nifti voxdim
#' @importFrom stats median
#' @importFrom ANTsRCore resampleImage iMath smoothImage as.antsImage
#' @importFrom ANTsRCore antsImageClone
segment_lung = function(img,
                        mask = TRUE,
                        lthresh = -300,
                        verbose = TRUE

                        ) {

  if (verbose) {
    message("Checking Inputs")
  }
  reg_img = check_ants(img)
  img = antsImageClone(reg_img)
  vres = voxdim(reg_img)
  if (verbose) {
    message("Resampling Image")
  }
  reg_img = resampleImage(reg_img, c(1,1,1))


  ##############################
  # 1024 should be lower limit
  ##############################
  if (verbose) {
    message("Making Positive Values")
  }
  adder = 1025
  reg_img = as.array(reg_img)

  reg_img = reg_img + adder
  if (verbose) {
    message("Setting voxel ranges")
  }
  reg_img[reg_img < 0] = 0
  reg_img[reg_img > 3071 + adder] = 3071 + adder
  reg_img = as.antsImage(reg_img, reference = img)

  if (verbose) {
    message("# Getting Humans")
  }
  ss = iMath(reg_img, "PeronaMalik", 10, 5)
  body = ss > (0 + adder)
  body = iMath(img = body, operation = "GetLargestComponent")
  inds = getEmptyImageDimensions(body)
  ss = maskEmptyImageDimensions(
    img = ss,
    inds = inds,
    mask.value = 0)


  body = ss > (lthresh + adder)
  body = iMath(img = body, operation = "GetLargestComponent")
  ebody = coarse_body(body)
  # inds = getEmptyImageDimensions(body)
  # ss = maskEmptyImageDimensions(
  #   img = ss,
  #   inds = inds,
  #   mask.value = 0)
  ss = mask_img(ss, ebody)

  # if (verbose) {
  #   message("# smoothing image")
  # }
  # ss = smoothImage(inimg = ss, sigma = 10,
  #                  max_kernel_width = 200)
  # # smoothed image is greater than some thresh
  # if (mask) {
  #   vals = as.numeric(ss)
  # } else {
  #   vals = as.numeric(ss)
  #   vals = vals[ as.numeric(body) > 0]
  # }
  # med = median(vals)
  # body = ss > med
  #
  # # body = ss > (-100 + adder)
  # if (verbose) {
  #   message("# Getting Humans")
  # }
  #
  # # zero padding in case for any connectedness
  # zp = as.array(body)
  # kdim = c(1,1,1)
  # zp = zero_pad(zp, kdim = kdim)
  # zp = as.antsImage(zp, reference = body)
  # # getting connected component
  # cc = iMath(img = zp, operation = "GetLargestComponent")
  # rm(list = "zp"); gc(); gc()
  # cc = as.array(cc)
  # cc = zero_pad(cc, kdim = kdim, invert = TRUE)
  # cc = as.antsImage(cc, reference = body)
  #
  # if (verbose) {
  #   message("# Filling Holes")
  # }
  # # cc = iMath(img = cc, operation = "FillHoles")
  # cc = filler(cc, fill_size = 60)
  # # cc = filler(cc, fill_size = 40)
  # cc = filler(cc, fill_size = 5, dilate = FALSE)
  #
  # # Dropping non-human stuff
  # inds = getEmptyImageDimensions(cc)
  #
  # if (verbose) {
  #   message("# Making New image")
  # }
  #
  # # Don't want to drop the indices,
  # # just blank them out
  # newimg = maskEmptyImageDimensions(img = reg_img, inds = inds)
  newimg = ss
  cc = ebody

  lung = newimg < (lthresh + adder) & newimg > 0 & cc == 1
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
  lung = mask_img(img, lung_mask)

  # reg_img = reg_img - adder
  L = list(img = img,
           lung_mask = lung_mask,
           lung = lung
  )
  return(L)
}
