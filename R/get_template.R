#' Create a lung template from warped images and composite transformations
#'
#' @param folder_warp Folder path for warped images
#' @param folder_comp Folder path for composite transformations
#' @param outfile File path to save new template
#' @param gradientStep Gradient step size
#' @param mask Logical statement. TRUE if outfile should be binary.
#' @param verbose Print output messages
#'
#' @return Template and its filename
#' @importFrom ANTsRCore antsAverageImages antsImageWrite antsImageRead antsApplyTransforms iMath
#' @export
#'
#' @examples
get_template = function(folder_warp,
                        folder_comp,
                        outfile,
                        gradientStep = 0.2,
                        mask = TRUE,
                        verbose = TRUE)
{
  if(verbose){
    message("Averaging warped images")
  }
  lf = list.files(folder_warp)
  template = antsAverageImages(paste0(folder_warp,lf))

  if(verbose){
    message("Averaging composite transformations")
  }
  lf = list.files(folder_comp)
  avg_comp = antsAverageImages(paste0(folder_comp,lf)) * (-1 * gradientStep)
  fn_temp = tempfile(fileext = ".nii.gz")
  antsImageWrite(avg_comp, fn_temp)

  if(verbose){
    message("Applying average transformation to average image")
  }
  template = antsApplyTransforms(template, template, fn_temp)

  if(verbose){
    message("Smoothing")
  }
  template = template * 0.5 + iMath(template, "Sharpen") * 0.5
  if(mask){
    template = template >= 0.5
  }

  if(verbose){
    message("Saving new template")
  }
  antsImageWrite(template, outfile)
  return(list(transform = template, transform_file = outfile))
}