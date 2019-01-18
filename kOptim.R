
#----------------------------------------------------------------------------------
# K OPTIM
#----------------------------------------------------------------------------------
#
# name:         kOptim.R
# description:  search optim k to GGC measure using linear regressions  
# input:  - time series x and y   
#         - k value that is the maximun lag in the AR(p) models
#         - kInf value, is the is the lower limit calculated for the lag k 
# output:	max value of k using BIC criteria, in the regressions: 
#                     x_0 ~ x1 +...+ xk + y0 +...+ yk
#                     y_0 ~ x0 +...+ xk + y1 +...+ yk

# 
# creadet by:   Carolina Gamboa

#----------------------------------------------------------------------------------




kOptim <- function(x, y, k, kinf){
  
  
  N <- length(x) # size time series 
  M_d <- matrix(nrow = N-k, ncol = 2*k+2)
  
  for(i in 1:(k+1)){
    M_d[, i]       <- x[i:(N-k+i-1)]
    M_d[, (i + k+ 1)] <- y[i:(N-k+i-1)]
  } 
  
  M_d <- data.frame( M_d)
  names(M_d) <- c(paste("x.lag", k:0, sep = ""), 
                  paste("y.lag", k:0, sep = ""))

  # regression : x_0 ~ x_1 +...+ x_jj + y_0 +...+ y_jj, and selection of the maximun order (p1) 
  # using BIC criteria 
  # the jj value is between 0 and k value. (see inputs in the header)
  
  bic1 <- data.frame()

  for (jj in kinf:(k)){
    b1 <- BIC(model <- lm(x.lag0 ~ . , data = M_d[, c((k+1-jj):(k+1), 
                                                      (2*(k+1)-jj):(2*(k+1)))]))
    bic1 <- rbind(bic1, c(jj, b1))
  }
  
  names(bic1) <- c("k", "BIC")
  p1 <- bic1[which.min(bic1$BIC), "k"]
  
  # if in the linear regressions the maximun order p1 is equal to k we stop the function.  
  # otherwise (p1 < k) the next part determines the maximum order using the regressions
  # y_0 ~ y_1 +...+ y_jj + x_0 +...+ x_jj, 
  # and selection of the order p2 using BIC criteria
  
  
  if(p1 < k){
    bic2 <- data.frame()
    
    for (jj in kinf:(k)){
      b2 <- BIC(model <- lm(y.lag0 ~ . , data = M_d[, c((k+1-jj):(k+1), 
                                                        (2*(k+1)-jj):(2*(k+1)))]))
      bic2 <- rbind(bic2, c(jj, b2))
    }
    
    names(bic2) <- c("k", "BIC")
    p2 <- bic1[which.min(bic2$BIC), "k"]
    p  <- max(p1, p2)
    
  }
  
  if(!exists("p")) p <- p1
  
  return(p)
}

