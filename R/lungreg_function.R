reg_helper = function(
  moving_mask,
  fixed_mask,
  value = 1,
  moving = NULL,
  verbose = FALSE,
  outprefix = NULL,
  add_prefix = "_left",
  mask_interpolator = "nearestNeighbor",
  interpolator = "linear",
  save_memory = TRUE,
  ...
) {

  if (is.null(outprefix)) {
    outprefix = tempfile()
  }

  # Just the left lobe
  # moving_mask = moving_mask == value | moving_mask == 3
  # fixed_mask = fixed_mask == value  | fixed_mask == 3
  moving_mask = moving_mask == value
  fixed_mask = fixed_mask == value

  # if any voxels are in left mask
  if (sum(moving_mask) > 0) {
    reg = antsRegistration(
      fixed = fixed_mask,
      moving = moving_mask,
      typeofTransform = "SyN",
      verbose = verbose,
      ...)

    #compose the transforms
    outfile = paste0(outprefix, add_prefix, ".nii.gz")
    composed = antsApplyTransforms(
      fixed = fixed_mask,
      moving = moving_mask,
      transformlist = reg$fwdtransforms,
      compose = tempfile())
    file.copy(composed, outfile, overwrite = TRUE)

    transformed_mask = antsApplyTransforms(
      fixed = fixed_mask,
      moving = moving_mask,
      transformlist = outfile,
      interpolator = mask_interpolator)

    reg$warpedmovout = NULL
    reg$warpedfixout = NULL
    for (i in 1:10) {
      gc();
    }

    if (!is.null(moving)) {
      moving = check_ants(moving)
      image = moving * moving_mask
      transformed = antsApplyTransforms(
        fixed = fixed_mask,
        moving = image,
        transformlist = outfile,
        interpolator = interpolator)
    } else {
      transformed = NULL
      image = NULL
    }
  } else {
    return(NULL)
  }
  L = list(composed_fwdtransforms = outfile)
  L$transformed = transformed
  L$transformed_mask = transformed_mask
  if (!save_memory) {
    L$masked_image = image
    L$moving_mask = moving_mask
    L$fixed_mask = fixed_mask
  } else {
    rm(list = c("image", "moving_mask", "fixed_mask"))
    for (i in 1:10) {
      gc();
    }
  }
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
#' @param sides Do both left and right or only one?
#' @param outprefix Path to put the transformations,
#' passed to \code{\link{antsRegistration}}
#' @param verbose Print diagnostic messages
#' @param interpolator interpolator used to apply transformation to image,
#' passed to \code{\link{antsApplyTransforms}}
#' @param mask_interpolator interpolator used to apply transformation to moving mask,
#' passed to \code{\link{antsApplyTransforms}}
#' @param save_memory do garbage collection and not save intermediate
#' files in registration.
#' @param ... addition arguments to pass to \code{\link{antsRegistration}}
#'
#' @return A list of registrations, transformed images, each for
#' left and right separately
#' @importFrom ANTsRCore antsRegistration antsApplyTransforms
#' @export
register_lung_mask = function(
  moving_mask,
  fixed_mask,
  moving = NULL,
  sides = c("left", "right"),
  outprefix = NULL,
  verbose = FALSE,
  mask_interpolator = "nearestNeighbor",
  interpolator = "linear",
  save_memory = TRUE,
  ...
) {

  sides = match.arg(sides, several.ok = TRUE)
  # verbose = FALSE;
  moving_mask = check_ants(moving_mask)
  fixed_mask = check_ants(fixed_mask)


  if (is.null(outprefix)) {
    outprefix = tempfile()
  }
  if (!is.null(moving)) {
    moving = check_ants(moving)
  }


  args = list(
    moving_mask = moving_mask,
    fixed_mask = fixed_mask,
    moving = moving,
    verbose = verbose,
    outprefix = outprefix,
    mask_interpolator = mask_interpolator,
    interpolator = interpolator,
    save_memory = save_memory,
    ...)

  if ("left" %in% sides) {
    # run left
    args$add_prefix = "_left"
    args$value = 1
    reg_left = do.call("reg_helper", args = args)
  } else {
    reg_left = NULL
  }

  if ("right" %in% sides) {
    # run right
    args$add_prefix = "_right"
    args$value = 2
    reg_right = do.call("reg_helper", args = args)
  } else {
    reg_right = NULL
  }


  res = list()
  res$reg_left = reg_left
  res$reg_right = reg_right
  res$interpolator = interpolator
  res$mask_interpolator = mask_interpolator
  res$outprefix = outprefix
  return(res)
}