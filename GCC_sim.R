
#----------------------------------------------------------------------------------
# GCC_d
#----------------------------------------------------------------------------------

# name:         GCC_sim.R
# description:  calcule the values of GCC measure
# inputs:  
# - two time series
# - lag k , if the user don't know this value, we use by default floor(N/10).
# output:       values of GCC similarity measure
#
# creadet by:   Carolina Gamboa



GCC_sim <- function(x, y, k){
  N <- length(x)
  
  if(missing(k)){
    k <- floor(N/10)
  }

  M_xy <- matrix(nrow = N-k, ncol = 2*(k+1))
  
  for(i in 1:(k+1)){
    M_xy[, i]       <- x[i:(N-k+i-1)]
    M_xy[, i+(k+1)] <- y[i:(N-k+i-1)]
  }
  
  
  M_x <- M_xy[, 1:(k+1)]
  M_y <- M_xy[, (k+2):(2*(k+1))]
  
  # Sample corr
  
  R_xy <- cor(M_xy)
  
  if(is.null(dim(M_x))){
    R_x  <- cor(M_x, M_x)
    R_y  <- cor(M_y, M_y)
    GCC <- 1 - det(R_xy)^(1/ (1 * (k+1)))  / (R_x^(1/ (1 * (k+1))) * R_y^(1/ (1 * (k+1))))
    
  } else  { 
    R_x  <- cor(M_x)
    R_y  <- cor(M_y)
    GCC <- 1 - det(R_xy)^(1/ (1 * (k+1)))  / (det(R_x)^(1/ (1 * (k+1))) * det(R_y)^(1/ (1 * (k+1))))
    
  }
  
  
  return(GCC)
  
}
