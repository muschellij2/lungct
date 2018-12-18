#' Register lungs, with options to save warped image and
#' composed transformation
#'
#' @param infile_template File path for the template (fixed) image
#' @param infile_moving File path for the moving image
#' @param template_masked Logical. Is the template already masked?
#' @param mask_value Value to mask moving image, and template if necessary
#' @param typeofTransform Type of transformation for registration,
#' see \code{\link{antsRegistration}} for options
#' @param outfile_warp File name to save warped image. If NULL, warped image is not saved to file.
#' @param outfile_comp File name to save composite transformation. If
#' \code{NULL}, transformation is not saved to file.
#' @param interpolator interpolation done after applying transformation,
#' passed to \code{\link{antsApplyTransforms}}
#' @param verbose Print messages.  If \code{> 1}, then more verbosity is
#' output.
#' @param ... additional arguments to pass to \code{\link{antsRegistration}}
#'
#' @return A list with the registration output and warps.
#'
#' @importFrom extrantsr check_ants
#' @importFrom ANTsRCore antsImageRead antsRegistration antsImageWrite antsImageClone antsApplyTransforms
#' @export
register_lung = function(
  infile_template,
  infile_moving,
  template_masked = TRUE,
  mask_value = c(1, 2),
  typeofTransform = "SyN",
  outfile_warp = NULL,
  outfile_comp = NULL,
  verbose = TRUE,
  interpolator = "nearestNeighbor",
  ...) {
  if (verbose) {
    message("# Checking/Reading Imaging data")
  }
  orig_fixed = check_ants(infile_template)
  orig_moving = check_ants(infile_template)
  if (verbose) {
    message("# Creating directories for outfiles if necessary")
  }
  if (!is.null(outfile_warp)) {
    dir.create(dirname(outfile_warp),
               showWarnings = FALSE,
               recursive = TRUE)
  }
  if (!is.null(outfile_comp)) {
    dir.create(dirname(outfile_comp),
               showWarnings = FALSE,
               recursive = TRUE)
  }

  if (verbose) {
    message("# Masking images")
  }

  ####################################
  # Allow this to run for both left and right
  ####################################
  L = vector(mode = "list", length = length(mask_value))
  for (irun in seq_along(mask_value)) {
    imask_value = mask_value[irun]
    if (!template_masked) {
      fixed = orig_fixed == mask_value
    } else {
      fixed = antsImageClone(orig_fixed)
    }
    moving = orig_moving == mask_value

    if (verbose) {
      message("# Registering images")
    }
    reg = antsRegistration(
      fixed = fixed,
      moving = moving,
      typeofTransform = typeofTransform,
      verbose = verbose > 1,
      ...
    )

    if (!is.null(outfile_warp)) {
      if (verbose) {
        message("# Saving warped image to file")
      }
      antsImageWrite(reg$warpedmovout, outfile_warp)
    }

    if (!is.null(outfile_comp)) {
      if (verbose) {
        message("# Changing pixel type to double")
      }
      m1 = antsImageClone(moving, out_pixeltype = "double")
      f1 = antsImageClone(fixed, out_pixeltype = "double")

      if (verbose) {
        message("# Composing transformations")
      }
      fn_temp = antsApplyTransforms(
        fixed = f1,
        moving = m1,
        transformlist = reg$fwdtransform,
        compose = reg$fwdtransform[1],
        interpolator = interpolator,
        verbose = verbose > 1
      )
      composed = antsImageRead(fn_temp)

      if (verbose) {
        message("# Saving composed transformation to file")
      }
      antsImageWrite(composed, outfile_comp)
      rm(fn_temp)
      rm(m1)
      rm(f1)
      rm(composed)
    }
    sub_L = list(registration = reg)
    sub_L$warp = outfile_warp
    sub_L$composed = outfile_comp
    sub_L$mask_value = imask_value
    L[[irun]] = sub_L
    rm(fixed)
    rm(moving)
    rm(reg)
    for (i in 1:10) gc()
  }
  if (length(mask_value) == 1) {
    L = L[[1]]
  }
  return(L)
}