template_helper = function(folder_warp,
                           folder_comp,
                           sides = sides,
                           gradientStep = 0.2,
                           mask = TRUE,
                           verbose = TRUE)
{
  if(verbose){
    message(paste0("Averaging warped images for ",sides," lung"))
  }
  lf = list.files(folder_warp)
  template = antsAverageImages(paste0(folder_warp,lf))

  if(verbose){
    message(paste0("Averaging composite transformations for ",sides," lung"))
  }
  lf = list.files(folder_comp)
  avg_comp = antsAverageImages(paste0(folder_comp,lf)) * (-1 * gradientStep)
  fn_temp = tempfile(fileext = ".nii.gz")
  antsImageWrite(avg_comp, fn_temp)

  if(verbose){
    message(paste0("Applying average transformation to average image for ",sides," lung"))
  }
  template = antsApplyTransforms(template, template, fn_temp)

  if(verbose){
    message(paste0("Smoothing ",sides," lung"))
  }
  template = template * 0.5 + iMath(template, "Sharpen") * 0.5
  if (mask){
    template = template >= 0.5
  }

  return(template)
}

#' Lung Template Creation
#'
#' Create a new lung template from warped images and composite transformations. If the DSC < 0.99, more iterations should be performed.
#'
#' @param folder_warp Folder path for warped images
#' @param folder_comp Folder path for composite transformations
#' @param sides Do both left and right or only one?
#' @param gradientStep Gradient step size
#' @param mask Logical statement. TRUE if template should be binary.
#' @param verbose Print output messages
#'
#' @return New Template. Right lung = 1, left lung = 2, non-lung = 0.
#' @importFrom ANTsR antsAverageImages antsImageWrite antsImageRead iMath antsApplyTransforms
#' @export
get_template = function(folder_warp,
                        folder_comp,
                        sides = c("right", "left"),
                        gradientStep = 0.2,
                        mask = TRUE,
                        verbose = TRUE)
{

  sides = match.arg(sides, several.ok = TRUE)
  if (!file.exists(folder_warp)) {
    stop("Folder path for warped images does not exist")
  }
  if (!file.exists(folder_comp)) {
    stop("Folder path for transformations does not exist")
  }

  args = list(
    folder_warp,
    folder_comp,
    gradientStep = gradientStep,
    mask = mask,
    verbose = verbose)

  if ("right" %in% sides) {
    # run right
    args$sides = "right"
    template_right = do.call("template_helper", args = args)
  } else {
    template_right = NULL
  }

  if ("left" %in% sides) {
    # run left
    args$sides = "left"
    template_left = do.call("template_helper", args = args)
  } else {
    template_left = NULL
  }

  res = list()
  res$template_right = template_right
  res$template_left = template_left
  res$template = template_right + 2*template_left
  return(res)

}