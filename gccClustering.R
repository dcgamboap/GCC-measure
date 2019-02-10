

#----------------------------------------------------------------------------------
# gccClustering
#----------------------------------------------------------------------------------
#
# name:         gccClustering.R
# description:  In this program a hierarchical grouping is constructed, 
#               increasing the efficiency of the gap and silhouette statistics, 
#               cleaning the database of outliers.This implementation was taken 
#               from the developments made by Professor Andres Alonso.  

# input: 
        # - zData: time series
        # - zLag : selected k  
        # - Percentage: % of observations that are considered a small group.
        # - threshold: percentage of observations that are not considered outliers. 
        # - toPlot: receives TRUE or FALSE values. If the value is TRUE, 
        #           the dendograms will be printed before and after "cleaning" the outliers.
#                
# output:	
          # - R: assignments of time series to each of the groups. 
          # - NCL: number of groups

# 
# creadet by:   Andrés Alonso, modified by Diana Gamboa. 
#----------------------------------------------------------------------------------


#----------------------------------------------------------------------------------
# function and library
#----------------------------------------------------------------------------------
# library

#setwd("C:/Users/Carolina/Documents/2018/TESIS/programa final/programas/function")

library('forecast')
library('tseries')
library('cluster')

source("GCCmatrix.R")
source("silhouetteNClus.R")
source("gap.R")

#----------------------------------------------------------------------------------
# main
#----------------------------------------------------------------------------------
gccClustering <- function(zData, zLag, Percentage, Threshold, toPlot){
  
  if(missing(zLag)) DM <- GCCmatrix(zData)
  if(missing(Percentage)) Percentage <- 0.05
  if(missing(Threshold)) Threshold <- 0.9
  if(missing(toPlot)) toPlot <- FALSE  
  
  PP <- ncol(zData)
  DM  <- GCCmatrix(zData, zLag)
  
  zLag <- DM$k
  distanceMatrix <- as.dist(DM$DM)
  
  lr <- hclust(distanceMatrix, method = "single")

  PRC <- quantile(lr$height, probs = Threshold)
  where  <- which(lr$height>PRC )
  
  if(length(where)>0){
    ncMax <- max(PP - where[1], 20)
  } else {
    ncMax <- 20
  }
  
  clusters <- cutree(lr, k = ncMax)
  TAB <- table(clusters)
  
  W1 <- which(TAB <= Percentage*PP)
  W2 <- which(TAB > Percentage*PP)
  
  R1 <- NULL
  R2 <- NULL
  for(i in 1:length(W1)){
    R1 <- c(R1, which(clusters == W1[i]))
  }
  for(i in 1:length(W2)){
    R2 <- c(R2, which(clusters == W2[i]))
  }
  
  if(toPlot==TRUE){
    tab <- data.frame(labels = names(TAB), 
                      abs.freq = as.matrix(TAB)[, 1], 
                      rel.freq = as.matrix(TAB)[, 1]/sum(TAB))
    rownames(tab) <- NULL
    cat("> first cluster assignation \n\n")
    print(tab)
    cat("\n")
    
    par(mfrow = c(2, 1))
    plot(lr, hang = -1)
  }
  
  
  D <- DM$DM
  
  D <- D[R2,R2]
  
  lr <- hclust(as.dist(D), "single")

  ncMax <- 20
  clusters <- cutree(lr, 1:ncMax)
  
  NCL1 <- silhouetteNClus(ncMax, D, "single")
  NCL2 <- GAPdistance(D, clusters, 100)
  
  if(NCL1 != NCL2)
    cat("[Silhouette = ", NCL1, "GAP = ", NCL2, "]", "\n\n")
  
  NCL <- NCL1
  
  
  TAB <- table(clusters[, NCL]) 
  
  
  
  W1 <- which(TAB <= Percentage*PP)
  W2 <- which(TAB > Percentage*PP)
  
  R12 <- NULL
  R22 <- NULL
  for(i in 1:length(W1)){
    R12 <- c(R12, which(clusters[, NCL] == W1[i]))
  }
  for(i in 1:length(W2)){
    R22 <- c(R22, which(clusters[, NCL] == W2[i]))
  }
  
  CL <- rep(0, PP)
  CL[R2[R22]] <- clusters[R22, NCL]
  
  where  <- which(CL == 0)
  
  
  nearest = 1
  while(length(where)>0){
    nearest = nearest +1
    for(i in 1:length(where)){
      I1 <- order(DM$DM[where[i], ])
      CL[where[i]] <- CL[I1[nearest]]
    }
    where <- which(CL == 0)
  }
  
  TAB <- table(CL)
  
  if(toPlot==TRUE){
    tab <- data.frame(labels = as.numeric(as.character(names(TAB))), 
                      abs.freq = as.matrix(TAB)[, 1], 
                      rel.freq = as.matrix(TAB)[, 1]/sum(TAB))
    rownames(tab) <- NULL
    cat("> final cluster assignation \n \n")
    print(tab)
    cat("\n\n")
    
    plot(lr, hang = -1)
  }
  
  NCL <- sum(as.numeric(names(TAB))>0) - sum(TAB==0)
  if(max(as.numeric(names(TAB)))>NCL){
    WCL <- CL
    k <- 1
    for(i in 1:max(as.numeric(names(TAB)))){
      if(!is.na(tab[i, 1])){ 
        C <- as.numeric(names(TAB))[i]
        WCL[CL==C] <- k
        k <- k + 1
        }
      }
    
    CL <- WCL
  }
  
    
  return(CL)

}






