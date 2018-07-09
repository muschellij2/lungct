reg_helper = function(moving_mask,
                      fixed_mask, value = 1,
                      moving = NULL,
                      verbose = FALSE,
                      outprefix = NULL,
                      add_prefix = "_left_",
                      interpolator = "linear") {

  if (is.null(outprefix)) {
    outprefix = tempfile()
  }

  # Just the left lobe
  moving_mask = moving_mask == value | moving_mask == 3
  fixed_mask = fixed_mask == value  | fixed_mask == 3

  # if any voxels are in left mask
  if (sum(moving_mask) > 0) {
    reg = antsRegistration(
      fixed = fixed_mask,
      moving = moving_mask,
      typeofTransform = "SyN",
      verbose = verbose,
      outprefix = paste0(outprefix, add_prefix))
    if (!is.null(moving)) {
      moving = check_ants(moving)
      image = moving * moving_mask
      transformed= antsApplyTransforms(
        fixed = fixed_mask,
        moving = image,
        transformlist = reg$fwdtransforms,
        interpolator = interpolator)
    } else {
      transformed = NULL
      image = NULL
    }
  } else {
    return(NULL)
  }
  L = list(reg_out = reg)
  L$transformed = transformed
  L$masked_image = image
  L$moving_mask = moving_mask
  L$fixed_mask = fixed_mask
  return(L)
}

#' Lung Registration for Left and Right Lungs
#'
#' @param moving_mask Mask of moving image, left has values
#' 1 and 3 (3 if both left and right), right has values 2 and 3
#' @param fixed_mask Mask of fixed image, left has values
#' 1 and 3 (3 if both left and right), right has values 2 and 3
#' @param moving image of moving image, usually CT, will be masked
#' by \code{moving_mask}
#' @param outprefix Path to put the transformations,
#' passed to \code{\link{antsRegistration}}
#' @param verbose Print diagnostic messages
#' @param interpolator interpolator used to apply transformation to image,
#' passed to \code{\link{antsApplyTransforms}}
#'
#' @return A list of registrations, transformed images, each for
#' left and right separately
#' @importFrom ANTsRCore antsRegistration antsApplyTransforms
#' @export
register_lung_mask = function(moving_mask,
                              fixed_mask,
                              moving = NULL,
                              outprefix = NULL,
                              verbose = FALSE,
                              interpolator = "linear"
) {
  # verbose = FALSE;
  moving_mask = check_ants(moving_mask)
  fixed_mask = check_ants(fixed_mask)


  if (is.null(outprefix)) {
    outprefix = tempfile()
  }
  if (!is.null(moving)) {
    moving = check_ants(moving)
  }

  reg_left = reg_helper(
    moving_mask,
    fixed_mask,
    value = 1,
    moving = moving,
    verbose = verbose,
    outprefix = outprefix,
    add_prefix = "_left_",
    interpolator = interpolator)

  reg_right = reg_helper(
    moving_mask,
    fixed_mask,
    value = 2,
    moving = moving,
    verbose = verbose,
    outprefix = outprefix,
    add_prefix = "_right_",
    interpolator = interpolator)
  res = list()
  res$reg_left = reg_left
  res$reg_right = reg_right
  res$interpolator = interpolator
  res$outprefix = outprefix
  return(res)
}