#' Register lungs, with options to save warped image and composed transformation
#'
#' @param infile_template File path for the template (fixed) image
#' @param infile_moving File path for the moving image
#' @param template_masked Logical. Is the template already masked?
#' @param mask_value Value to mask moving image, and template if necessary
#' @param typeofTransform Type of transformation for registration
#' @param outfile_warp File name to save warped image. If NULL, warped image is not saved to file.
#' @param outfile_comp File name to save composite transformation. If NULL, transformation is not saved to file.
#' @param out Logical. Should the registration object be returned?
#' @param verbose Print messages
#' @param verbose_reg Print messages, specifically related to registration function
#'
#' @return antsRegistration object, if out = TRUE
#' @importFrom ANTsRCore antsImageRead antsRegistration antsImageWrite antsImageClone antsApplyTransforms
#' @export
register_lung = function(infile_template,
                         infile_moving,
                         template_masked = TRUE,
                         mask_value = 1,
                         typeofTransform = "SyN",
                         outfile_warp = NULL,
                         outfile_comp = NULL,
                         out = FALSE,
                         verbose = TRUE,
                         verbose_reg = FALSE)
{

  if(verbose){
    message("# Checking if folder/file paths exist")
  }
  try(!file.exists(infile_template), stop("File path for template image does not exist"))
  try(!file.exists(infile_moving), stop("File path for moving image does not exist"))
  if(!is.null(outfile_warp)){
    try(!file.exists(dirname(outfile_warp)), stop("Folder path for new warped image does not exist"))
  }
  if(!is.null(outfile_comp)){
    try(!file.exists(dirname(outfile_comp)), stop("Folder path for new composite transformation does not exist"))
  }


  if(verbose){
    message("# Reading in images")
  }
  fixed = antsImageRead(infile_template)
  moving = antsImageRead(infile_moving)

  if(verbose){
    message("# Masking images")
  }
  if(!template_masked){
    fixed = fixed == mask_value
  }
  moving = moving == mask_value

  if(verbose){
    message("# Registering images")
  }
  reg = antsRegistration(fixed = fixed,
                         moving = moving,
                         typeofTransform = typeofTransform,
                         verbose = verbose_reg)

  if(!is.null(outfile_warp)){
    if(verbose){
      message("# Saving warped image to file")
    }
    antsImageWrite(reg$warpedmovout, outfile_warp)
  }

  if(!is.null(outfile_comp)){

    if(verbose){
      message("# Changing pixel type to double")
    }
    m1 = antsImageClone(moving, out_pixeltype = "double")
    f1 = antsImageClone(fixed, out_pixeltype = "double")


    if(verbose){
      message("# Composing transformations")
    }
    fn_temp = antsApplyTransforms(fixed = f1,
                                  moving = m1,
                                  transformlist = reg$fwdtransform,
                                  compose = reg$fwdtransform[1])
    composed = antsImageRead(fn_temp)

    if(verbose){
      message("# Saving composed transformation to file")
    }
    antsImageWrite(composed, outfile_comp)
  }

  if(out){
    return(reg)
  } else{
    return(NULL)
  }

}