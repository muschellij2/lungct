reg_helper = function(
  moving_mask,
  fixed_mask,
  value = 1,
  moving = NULL,
  composeTransforms = NULL,
  add_prefix = "_right_",
  verbose = FALSE,
  typeofTransform = "SyN",
  mask_interpolator = "nearestNeighbor",
  interpolator = "linear",
  ...
) {


  # Which lobe
  moving_mask = moving_mask == value
  fixed_mask = fixed_mask == value

  # Register lung mask
  if (sum(moving_mask) > 0) {
    reg = antsRegistration(
      fixed = fixed_mask,
      moving = moving_mask,
      typeofTransform = typeofTransform,
      verbose = verbose,
      ...)

    # If composeTransforms in non-null
    if (!is.null(composeTransforms)){
      output_prefix = paste0(composeTransforms, add_prefix)
      fixed_mask = antsImageClone(fixed_mask, out_pixeltype = "double")
      moving_mask = antsImageClone(moving_mask, out_pixeltype = "double")
      composed = antsApplyTransforms(fixed = fixed_mask,
                                    moving = moving_mask,
                                    transformlist = reg$fwdtransforms,
                                    compose = output_prefix)
    } else {
      composed = NULL
    }

    # If moving image is non-null
    if (!is.null(moving)) {

      # Get images in correct format (otherwise antsApplyTransforms will throw errors)
      image = maskImage(moving, moving_mask)
      image = antsImageClone(image, out_pixeltype = "float")
      fixed_mask = antsImageClone(fixed_mask, out_pixeltype = "double")

      # Apply transformation to moving image
      transformed = antsApplyTransforms(
        fixed = fixed_mask,
        moving = image,
        transformlist = reg$fwdtransforms,
        interpolator = interpolator)

      # Mask with warped mask and template mask, just to be sure
      transformed = maskImage(transformed, reg$warpedmovout)
      transformed = maskImage(transformed, fixed_mask)
    } else {
      transformed = NULL
    }
  } else {
    return(NULL)
  }
  L = list(warped_mask = reg$warpedmovout,
           warped_img = transformed,
           fwdtransforms = reg$fwdtransforms,
           composedtransform = composed)
  return(L)
}

#' Lung Registration
#'
#' This function registers the right and left lung masks to a template mask. To register to the standard lung template mask, type \code{system.file("extdata", "lung_template_mask.nii.gz", package = "lungct").}
#'
#' @param moving_mask Mask of moving image. Right lung = 1, left lung = 2, non-lung = 0
#' @param fixed_mask Mask of fixed image. Right lung = 1, left lung = 2, non-lung = 0
#' @param moving Moving image to apply transformation
#' @param sides Choose to register right and/or left lungs.
#' @param verbose Print diagnostic messages
#' @param typeofTransform Type of transform, passed to \code{\link{antsRegistration}}
#' @param composeTransforms Prefix of output filename to save the composed forward transformations. The prefix will add comptx.nii.gz to the end.
#' @param mask_interpolator Interpolator used to apply transformation to moving mask,
#' passed to \code{\link{antsApplyTransforms}}
#' @param interpolator Interpolator used to apply transformation to image,
#' passed to \code{\link{antsApplyTransforms}}
#' @param ... addition arguments to pass to \code{\link{antsRegistration}}
#'
#' @return A list of warped masks, images, and transformations for
#' right and left lungs separately
#' @importFrom ANTsR antsRegistration antsApplyTransforms
#' @export
register_lung_mask = function(
  moving_mask,
  fixed_mask,
  moving = NULL,
  sides = c("right", "left"),
  verbose = FALSE,
  typeofTransform = "SyN",
  composeTransforms = NULL,
  mask_interpolator = "nearestNeighbor",
  interpolator = "linear",
  ...
) {

  sides = match.arg(sides, several.ok = TRUE)
  moving_mask = check_ants(moving_mask)
  fixed_mask = check_ants(fixed_mask)
  if (!is.null(moving)) {
    moving = check_ants(moving)
  }


  args = list(
    moving_mask = moving_mask,
    fixed_mask = fixed_mask,
    moving = moving,
    verbose = verbose,
    typeofTransform = typeofTransform,
    composeTransforms = composeTransforms,
    mask_interpolator = mask_interpolator,
    interpolator = interpolator,
    ...)


  if ("right" %in% sides) {
    # run right
    args$add_prefix = "_right_"
    args$value = 1
    reg_right = do.call("reg_helper", args = args)
  } else {
    reg_right = NULL
  }


  if ("left" %in% sides) {
    # run left
    args$add_prefix = "_left_"
    args$value = 2
    reg_left = do.call("reg_helper", args = args)
  } else {
    reg_left = NULL
  }


  res = list()
  res$right = reg_right
  res$left = reg_left
  return(res)
}
