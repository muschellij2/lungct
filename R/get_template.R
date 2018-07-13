template_helper = function(folder_warp,
                           folder_comp,
                           sides = sides,
                           outprefix = NULL,
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
  if(mask){
    template = template >= 0.5
    if(sides == "right"){
      template = template * 2
    }
  }

  if(verbose){
    message(paste0("Saving new template for ",sides," lung"))
  }
  outfile = paste0(outprefix,"_",sides,".nii.gz")
  antsImageWrite(template, outfile)
  return(list(template = template, template_file = outfile))
}

#' Create a lung template from warped images and composite transformations
#'
#' @param folder_warp Folder path for warped images
#' @param folder_comp Folder path for composite transformations
#' @param sides Do both left and right or only one?
#' @param outprefix Folder and file path to save new template (don't include extension)
#' @param gradientStep Gradient step size
#' @param mask Logical statement. TRUE if outfile should be binary.
#' @param verbose Print output messages
#'
#' @return Template and its filename
#' @importFrom ANTsRCore antsAverageImages antsImageWrite antsImageRead antsApplyTransforms iMath
#' @export
get_template = function(folder_warp,
                        folder_comp,
                        sides = c("left", "right"),
                        outprefix = NULL,
                        gradientStep = 0.2,
                        mask = TRUE,
                        verbose = TRUE)
{

  sides = match.arg(sides, several.ok = TRUE)
  try(!file.exists(folder_warp), stop("Folder path for warped images does not exist"))
  try(!file.exists(folder_comp), stop("Folder path for transformations does not exist"))
  try(!file.exists(dirname(outprefix)), stop("Folder path for new template does not exist"))

  args = list(
    folder_warp,
    folder_comp,
    outprefix = outprefix,
    gradientStep = gradientStep,
    mask = mask,
    verbose = verbose)

  if ("left" %in% sides) {
    # run left
    args$sides = "left"
    template_left = do.call("template_helper", args = args)
  } else {
    template_left = NULL
  }

  if ("right" %in% sides) {
    # run right
    args$sides = "right"
    template_right = do.call("template_helper", args = args)
  } else {
    template_right = NULL
  }

  res = list()
  res$template_left = template_left
  res$template_right = template_right
  return(res)

}