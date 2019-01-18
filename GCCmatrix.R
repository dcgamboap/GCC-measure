#----------------------------------------------------------------------------------
# clustering SV Model
#----------------------------------------------------------------------------------
#
# name:         GCCmatrix.R
# description:  built GCC similarity matrix based on dependences 
# input: 
# functions:  
            # serie_0 <- time series  matrix
            # k <- max dependence lag if the user knows this value. In other case the
            #        function calcules this value
            # model <- when the value k is unknown for the user, the specification of model
            #           determine if we works with GARCH model or AR model, if you don't specify 
            #             the model you will work with the default ARMA model.
#                
# output:	a list that contain: - a matrix object with the distance matrix
        #                      - the lag k used to calculate GCC measure
# creadet by:   Carolina Gamboa
#
#----------------------------------------------------------------------------------

#----------------------------------------------------------------------------------
# functions
#----------------------------------------------------------------------------------
# library

library('forecast')
library('tseries')

# this fuctions are used in the program. 
# the description of each function can be found in the header of the program. 
source("kOptim.R") 
source("GCC_sim.R")


#............................................................................
# Main  
#............................................................................


GCCmatrix <- function(serie_0, k, model){
  # first we determine the upper limit for "k" using AR(P) model
  
  if(missing(k)){
    
    if(missing(model)){
      
      setAR <- function(x){
        N <- length(x)
#        orderMax <- min(N-1, 10*log10(N))
        orderMax <- min(10, N/10)
        spfinal.bic <- Inf
        spfinal.order <- c(0,0,0)
        for(i in 0:orderMax){
          model0 <- arima(x, order=c(i, 0, 0))
          spcurrent.aic <- AIC(model0)
          spcurrent.bic <- BIC(model0)
          if (spcurrent.bic < spfinal.bic) {
            spfinal.bic <- spcurrent.bic
            spfinal.order <- c(i, 0, 0)
               }
          }
        arOrder <- spfinal.order[1]
        
        return(arOrder)
      }
      
      
      PP <- apply(serie_0, 2, setAR) # search max order in data set
      kSup <- max(PP) # upper limit for "k"
      
    } else{
      
      orderArch <- function(x){
        N <- length(x)
        #        orderMax <- min(N-1, 10*log10(N))
        orderMax <- min(5, N/10)
        aic <- NULL
        bic <- NULL
        for(i in 1:orderMax){
          model1 <- garchFit(substitute(~garch(i,0),list(i=i)), data = x)
          aic <- c(aic, (2*(i+2) + 2*model1@fit$value)/1000)
          bic <- c(bic, (log(1000)*(i+2) + 2*model1@fit$value)/1000)
        }
        return(which.min(bic))
      }
      
      # search max order in data set. 
      # Note that you are taking square root to the time series. 
      # I have only used this methodology when we want to calculate the GCC 
      # in the squares of the time series for heteroscedastic time series. 
      
      PP <- apply(sqrt(serie_0), 2, orderArch)  
      kSup <- max(PP) # upper limit for "k"
    }
    
    
    # linear regression
    # we search optim k using  linear regressions between 0 and kSup
    kOp <- data.frame()
    kOp1 <- 0
    kinf <- 0
    for(jj in 1:(ncol(serie_0)-1)){
      for (ii in (jj+1):ncol(serie_0)){
        if(kOp1 == kSup){
          break
        }
          kOp1 <- c(kOptim(serie_0[, jj], serie_0[, ii], kSup, kinf))
          
          kOp <- rbind(kOp, data.frame(i = jj, j = ii, kOp1))
          kinf <- kOp1
           # cat(jj, ii, ", ")
        
      }
    }
    
    
    if(nrow(kOp) != 0){ kDef <- max(kOp$kOp1)
    } else kDef= kSup
    
  } else kDef= k
  
  
  #------------------------------------------------------------------------------
  # similarity matrix
  
  nSerie <- ncol(serie_0)
  DM <- diag(x = 0, nrow = nSerie, ncol = nSerie)
  
  # construction of the dissimilarity matrix. 
  
  for(ii in 1:nSerie){
    for(jj in ii:nSerie){
      g <- GCC_sim(serie_0[, ii], serie_0[, jj], kDef)
      DM[ii, jj] <- 1 - g
      DM[jj, ii] <- 1 - g
    }
   # cat(ii, ",", " ")
  }
  
  colnames(DM) <- colnames(serie_0)
  
  sale <- list(DM = DM,   k_GCC = kDef)
  return(sale)
}

