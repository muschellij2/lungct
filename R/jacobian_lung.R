#' Calculate Jacobian Determinant Image for the Lungs
#'
#' This function calculates the Jacobian determinant image for the left and right lungs separately.
#'
#' @param right_transformation Transformations for the right lung
#' @param left_transformation Transformations for the left lung
#' @param mask Fixed mask (or standard lung template mask)
#' @param doLog Return the log Jacobian
#' @param geom Use the geometric Jacobian calculation
#' @param relative Returen the relative Jacobian (Jacobian divided by mean Jacobian)
#'
#' @return Jacobian Image
#' @importFrom ANTsR createJacobianDeterminantImage maskImage
jacobian_lung = function(
  right_transformation,
  left_transformation,
  mask,
  doLog = TRUE,
  geom = FALSE,
  relative = TRUE
) {

  # Create separate left and right masks
  right = mask == 1
  left = mask == 2

  # Create Jacobian Determinant Image
  jacob_right <- createJacobianDeterminantImage(domainImg = mask,
                                          tx = right_transformation,
                                          doLog = doLog,
                                          geom = FALSE)
  jacob_right <- maskImage(jacob_right, right)

  jacob_left <- createJacobianDeterminantImage(domainImg = mask,
                                                tx = left_transformation,
                                                doLog = doLog,
                                                geom = FALSE)
  jacob_left <- maskImage(jacob_left, left)

  # Put images together
  jacob = jacob_right + jacob_left

  if(relative){

    jacob_right_relative <- jacob_right/mean(jacob_right)
    jacob_left_relative <- jacob_left/mean(jacob_left)
    jacob_relative <- jacob_right_relative + jacob_left_relative

    mylist = list(jacob = jacob,
                  jacob_relative = jacob_relative)
  }else{
    mylist = jacob
  }

  return(mylist)
}