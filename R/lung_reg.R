#' @title CT Lung Registration
#' @description Wrapper of regsitration function for CT scans of the lung
#'
#' @param base Baseline scan - filename or \code{\link{"antsImage-class"}}
#' @param follow Follow up scan - filename or \code{\link{"antsImage-class"}}.
#' Registered to the baseline scan.
#' @param base_mask Baseline mask - binary image for baseline.
#' @param follow_mask Follow up mask - binary image for follow up scan
#' @param add_1025 Should 1025 be added to make values non-negative (recommended
#' for CT)
#' @param interpolator interpolation in baseline space,
#' passed to \code{\link{registration}}
#' @param typeofTransform transformation to use, passed to \code{\link{registration}}
#' @param ... additional options to pass to \code{\link{registration}}
#'
#' @return List of output, same as \code{\link{registration}}, with transforms
#' and output files, as well as an indicator for \code{add_1025} and
#' value added to the data
#' @export
#'
#' @importFrom extrantsr registration
lung_reg = function(base,
                    follow,
                    base_mask,
                    follow_mask,
                    add_1025 = TRUE,
                    interpolator = "LanczosWindowedSinc",
                    typeofTransform = "SyN",
                    ...) {

  L_base = reduce_scan(img = base, mask = base_mask)
  L_fup = reduce_scan(img = follow, mask = follow_mask)

  adder = function(L, add_val) {
    L$img = L$img + add_val
    L$img = maskImage(L$img, L$mask)
    return(L)
  }
  add_val = 0
  if (add_1025) {
    add_val = 1025
  }
  L_base = adder(L_base, add_val = add_val)
  L_fup = adder(L_fup, add_val = add_val)

  outprefix = tempfile()
  reg = registration(
    filename = L_fup$img,
    template.file = L_base$base,
    outprefix = outprefix,
    remove.warp = FALSE,
    interpolator = interpolator,
    typeofTransform = typeofTransform,
    ...
  )
  reg$add_1025 = add_1025
  reg$added_value = add_val
  return(reg)
}