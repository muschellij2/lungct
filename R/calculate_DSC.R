DSC <- function(x,y){
  x <- as.array(x)
  y <- as.array(y)
  tt <- sum( x &  y)
  tf <- sum( x & !y)
  ft <- sum(!x &  y)
  ff <- sum(!x & !y)
  tab = matrix(c(ff, tf, ft, tt), ncol=2)
  n = list(c("FALSE", "TRUE"), c("FALSE", "TRUE"))
  dimnames(tab) = n
  tab = as.table(tab)

  DSC = 2*tab[2,2] / (2*tab[2,2] + tab[1,2] + tab[2,1])
  return(DSC)
}
#' Dice Similarity Coefficient Calculation
#'
#' Calculate the Dice similarity coefficient (DSC) for two iterations of lung templates
#'
#' @param template0 Initial template mask. Right lung = 1, left lung = 2, non-lung = 0
#' @param template1 New template mask. Right lung = 1, left lung = 2, non-lung = 0
#' @param sides Calculate DSC on right and/or left lungs. Overall DSC is given by default.
#'
#' @return DSC on whole lung, right lung and left lung.
#' @export
calculate_DSC <- function(template0, template1, sides = c("right","left")){

  # Overall DSC
  x <- template0 > 0
  y <- template1 > 0
  DSC_overall <- DSC(x,y)


  if ("right" %in% sides){
    x <- template0 == 1
    y <- template1 == 1
    DSC_right <- DSC(x,y)
  }else{DSC_right = NULL}


  if ("left" %in% sides){
    x <- template0 == 2
    y <- template1 == 2
    DSC_left <- DSC(x,y)
  }else{DSC_left = NULL}

  return(list(DSC_overall = DSC_overall,
              DSC_right = DSC_right,
              DSC_left = DSC_left))
}
