register_lung = function(infile_template, 
                    infile_moving, 
                    template_masked = TRUE,
                    mask_value = 1, 
                    add_prefix = "_left",
                    typeofTransform = "SyN",
                    outfile_warp = NULL, 
                    outfile_comp = NULL,
                    out = FALSE,
                    verbose = TRUE)
{
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
                         typeofTransform = typeofTransform)
  
  if(!(outfile_warp == NULL)){
    if(verbose){
      message("# Saving warped image to file")
    }
    antsImageWrite(reg$warpedmovout, fn_warped)
  }

  if(!(outfile_comp == NULL)){
    
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
    antsImageWrite(composed, fn_composed)
  }
  
  if(out){
    return(reg)
  } else{
    return(NULL)
  }
  
}