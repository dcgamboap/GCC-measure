
#----------------------------------------------------------------------------------
# gccClustering
#----------------------------------------------------------------------------------
#
# name:         gccClustering.R

# description
# In this program a clustering of time series is made by using the procedure 
# in Alonso and Peña (2019). first the matrix of Generalized Cross correlation (GCC) is built by 
# using the subrutine GCCmatrix, then a hierarchical grouping is constructed and the number of 
#  clusters is selected by a modified silhouette statistics, as follows:
# (1) series that join the groups at a distance larger than the a given threshold of the 
#distribution of the distances are disregarded. 
# (2) A minimum  size for the groups is defined by Percentage, groups smaller than Percentage are 
# disregarded.
# (3) The final groups are obtained by applying the silhouette statistics to the part of the sample of 
# time series that verifies conditions (1) and (2).  
# (4) The series disregarded in steps (1) and (2) are asigned to the closest groups, assuming 
# that the series is not an outlier with respect the members of the group, as checked by the 
# median and the mad of the distances inside the group.

# input: 
#       - zData: matrix of time series in columns
#       - zLag : selected k  for the GCC (default is computed inside the program)
#       - threshold: percentail in the distribution of distances that define  
#           observations that are not considered outliers.
#       - Percentage: % relative size of the minimum group considered.
#       - toPlot: receives TRUE or FALSE values. If the value is TRUE, 
#           the dendograms will be printed before and after "cleaning" the outliers.
#     
#                
# output:	
#    -Tabla de number of clusters found and number of observations in each cluster.
#          Group 0 if it exists indicates the outlier group. 
#  - sal: is a list with two objects
#           - $labels: labels of clusters
#           - $groups: is a list of matrices. each matrix corresponds to the set of time 
#                       series that make up each group. 
#                       For example, $groups[[i]] contains the set of time series 
#                       that belong to the ith group. 
# - 

# 
# written by Andrés Alonso, Carolina Gamboa and Daniel Peña. 
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
gccClustering <- function(zData, zLag, Percentage, Threshold, toPlot, asignG){
  
  if(missing(zLag)){ 
    DM <- GCCmatrix(zData); cat("> k used for GCC:", DM$k,  "\n\n")
  }
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
        
  tab <- data.frame(labels = as.numeric(as.character(names(TAB))), 
                    abs.freq = as.matrix(TAB)[, 1], 
                    rel.freq = as.matrix(TAB)[, 1]/sum(TAB))
  rownames(tab) <- NULL        
  
  if(toPlot==TRUE){

 
    
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
  
  groupAsign <- list()
  for(i in 1:length(unique(CL))){
    groupAsign[[i]]  <-  zData[, which(CL == i)]
  }  
  

  cat("> number final of clusters", length(unique(CL)), "\n \n")
  
  cat("> frecuency table\n \n")

  TAB <- table(CL)
  tab <- data.frame(labels = as.numeric(as.character(names(TAB))), 
                    abs.freq = as.matrix(TAB)[, 1], 
                    rel.freq = as.matrix(TAB)[, 1]/sum(TAB))
 
  
  print(tab)
  cat("\n\n")
  
  sal <- list()
  sal$labels <- CL
  sal$groups <- groupAsign
  

  return(sal)

}
