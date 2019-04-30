
#----------------------------------------------------------------------------------
# MTSCluster
#----------------------------------------------------------------------------------
#
# name:         MTSCluster.R

# description
# In this program a clustering of time series is made by using the procedure 
# in Alonso and Peña (2019). first the matrix of Generalized Cross correlation (GCC) is built by 
# using the subrutine GCCmatrix, then a hierarchical grouping is constructed and the number of 
#  clusters is selected by a modified silhouette statistics, as follows:
# (1) series that join the groups at a distance larger than a given threshold of the 
# distribution of the distances are disregarded. 
# (2) A minimum  size for the groups is defined by Percentage, groups smaller than Percentage are 
# disregarded.
# (3) The final groups are obtained in two steps .First the silhouette statistics is applied to the set of  
# time series that verify conditions (1) and (2). Second, the series disregarded in steps (1) and (2) are candidates
# to be asigned to its closest group. It is checked using the median and the mad of the group if the point
# is or it is not an outlier with respect to the group. If it is an outlier it is included in a group 0 of outlier # series.
# The distance between a series and a group is usually to the closest in the group (simple linkage) but could be to
# the mean of the group

# input: 
#       - zData: matrix of time series in columns
#       - zLag : selected lag k for the GCC (default is computed inside the program)
#       - threshold: percentail in the distribution of distances that define  
#           observations that are not considered outliers.(default threshold value is 0.9)
#       - Percentage: % relative size of the minimum group considered.(The default value is  0.05.)
#       - toPlot: receives TRUE or FALSE values. If the value is TRUE, 
#           a clustermatrix plot of distances is presented.
#     
#                
# output:	
#        -Table of number of clusters found and number of observations in each cluster.
#              Group 0 if it exists indicates the outlier group. 
#        - sal: is a list with four objects
#           - $labels: assignments of time series to each of the groups. 
#           - $groups: is a list of matrices. each matrix corresponds to the set of time 
#                       series that make up each group. 
#                       For example, $groups[[i]] contains the set of time series 
#                       that belong to the ith group. 
#            $matrix    GCC distance matrix
#            $gmatrix    GCC distance matrices in each group 
#      A clustermatrix plot is included with the distances inside each group in the diagonal boxes and the distances #      between series in two groups in off-diagonal boxes
#  

# 
# written by Andrés Alonso, Carolina Gamboa and Daniel Peña. 
#----------------------------------------------------------------------------------



#----------------------------------------------------------------------------------
# function and library
#----------------------------------------------------------------------------------
# library

#setwd()  here the working directory should be defined 

library('forecast')
library('tseries')
library('cluster')

#source("GCCmatrix.R")
#source("silhouetteNClus.R")
#source("gap.R")
#source("graphMatrix.R")

#----------------------------------------------------------------------------------
# main
#----------------------------------------------------------------------------------
MTSCluster <- function(zData, zLag, Percentage, Threshold, toPlot){
  
  if(is.null(names(zData))){
    colnames(zData) <- 1:ncol(zData)
  }
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
    #plot(lr, hang = -1)
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
 
  TAB <- table(CL)
        
  tab <- data.frame(labels = as.numeric(as.character(names(TAB))), 
                    abs.freq = as.matrix(TAB)[, 1], 
                    rel.freq = as.matrix(TAB)[, 1]/sum(TAB))
  rownames(tab) <- NULL        
 
  
  NCL <- sum(as.numeric(names(TAB))>0) - sum(TAB==0)
  if(max(as.numeric(names(TAB)))>NCL){
    WCL <- CL
    k <- 1
    if(sum(names(TAB)==0)==1)     k <- 0

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
  for(i in 1:max(unique(CL))){
    groupAsign[[i]]  <-  zData[, which(CL == i)]
  }  
  
  groupAsign[["outlier"]] <- zData[, which(CL == 0)]

  
  groupDist <- list()
  for(i in 1:max(unique(CL))){
    groupDist[[i]]  <-  DM$DM[which(CL == i), which(CL == i)]
  }  
  
  groupDist[["outlier"]] <- DM$DM[which(CL == 0), which(CL == 0)]
  
  fMedian <- function(x) median(as.dist(x))
  fMAD <- function(x) mad(as.dist(x))
  
  medianG <- lapply(groupDist, fMedian)
  MADG <- lapply(groupDist, fMAD)

  criteriaG <- rep(0, max(CL))
  for(i in 1:max(CL)){
    criteriaG[i] <- medianG[[i]] + 3 * MADG[[i]]
  }
  
  names(criteriaG) <- 1:max(CL)
  
  # reasignation of outliers 
  
  where <- which(CL == 0)
  
  if(length(where) >0){
      nearbor <- DM$DM[where,-where]
      for(i in 1:nrow(nearbor)){
        disSort <- sort(nearbor[i, ])
        id <- names(disSort)[1]
        where0 <- which(colnames(zData)==id)
        groupId <- paste(CL[where0])
        if(disSort[1] < criteriaG[groupId]){
          CL[where[i]] <- CL[where0]
        }
      }

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
  sal$matrix <- DM$DM
  sal$gmatrix <- groupDist
  if(toPlot==TRUE){
    graphMatrix(DM$DM, CL)
  }

  return(sal)

}



#############################################################################################################
#
# name:         GCCmatrix.R
# description:  built the GCC similarity matrix between time series  
# input: 
# functions:  
# zData <- time series  matrix with series in columns
# k <- max dependence lag if the user knows this value. In other case the
#        function calcules this value
# model <- when the value k is unknown for the user, the specification of model
#           determine if we works with GARCH model or AR model, if you don't specify 
#             the model you will work with the default ARMA model.
#                
# output:	a list that contain: - a matrix object with the distance matrix
#                      - the lag k used to calculate GCC measure
# #
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


GCCmatrix <- function(zData, k, model){
  # first we determine the upper limit for "k" using AR(P) model
  
  if(missing(k)){
    
    if(missing(model)){
      # in setAR function we estimate AR models varying order p between 1 and orderMax parameter, 
      #and return the order that maximice the BIC criteria.  
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
      
      
      PP <- apply(zData, 2, setAR) # search max order in data set
      kSup <- max(PP) # upper limit for "k"
      
    } else{
      # in orderArch function we estimate ARCH(p) models varying order p between 1 and orderMax parameter, 
      #and return the order that maximice the BIC criteria.  
      
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
      
      PP <- apply(sqrt(zData), 2, orderArch)  
      kSup <- max(PP) # upper limit for "k"
    }
    
    
    # linear regression
    # we search optim k using  linear regressions between 0 and kSup
    kOp <- data.frame()
    kOp1 <- 0
    kinf <- 0
    for(jj in 1:(ncol(zData)-1)){
      for (ii in (jj+1):ncol(zData)){
        if(kOp1 == kSup){
          break
        }
        kOp1 <- c(kOptim(zData[, jj], zData[, ii], kSup, kinf))
        
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
  
  nSerie <- ncol(zData)
  DM <- diag(x = 0, nrow = nSerie, ncol = nSerie)
  
  # construction of the dissimilarity matrix. 
  
  for(ii in 1:nSerie){
    for(jj in ii:nSerie){
      g <- GCC_sim(zData[, ii], zData[, jj], kDef)
      DM[ii, jj] <- 1 - g
      DM[jj, ii] <- 1 - g
    }
    # cat(ii, ",", " ")
  }
  
  if(!is.null(colnames(zData))){
    colnames(DM) <- colnames(zData)    
  } else{
    colnames(DM) <- 1:ncol(zData)
    rownames(DM) <- 1:ncol(zData)
    
  }  
  sale <- list(DM = DM,   k_GCC = kDef)
  return(sale)
}


###############################################################################################################################################
#
# name:         GCC_sim.R
# description:  calcule the values of GCC measure
# inputs:  
# - two time series: xData and yData
# - lag k , if the user don't know this value, we use by default floor(N/10).
# output:       values of GCC similarity measure
#
#
###############################################################################################################################################

GCC_sim <- function(xData, yData, k){
  N <- length(xData)
  
  if(missing(k)){
    k <- floor(N/10)
  }
  
  M_xy <- matrix(nrow = N-k, ncol = 2*(k+1))
  
  for(i in 1:(k+1)){
    M_xy[, i]       <- xData[i:(N-k+i-1)]
    M_xy[, i+(k+1)] <- yData[i:(N-k+i-1)]
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

#----------------------------------------------------------------------------------
# K OPTIM
#----------------------------------------------------------------------------------
#
# name:         kOptim.R
# description:  search optim k to GGC measure using linear regressions  
# input:  - time series xData and yData   
#         - kMax value that is the maximun lag in the AR(p) models
#         - kInf value, is the is the lower limit calculated for the lag k 
# output:	max value of k using BIC criteria, in the regressions: 
#                     x_0 ~ x1 +...+ xk + y0 +...+ yk
#                     y_0 ~ x0 +...+ xk + y1 +...+ yk

#----------------------------------------------------------------------------------

kOptim <- function(xData, yData, kMax, kinf){
  
  
  N <- length(xData) # size time series 
  M_d <- matrix(nrow = N-kMax, ncol = 2*kMax+2)
  
  for(i in 1:(kMax+1)){
    M_d[, i]       <- xData[i:(N - kMax+i-1)]
    M_d[, (i + kMax + 1)] <- yData[i:(N - kMax+i-1)]
  } 
  
  M_d <- data.frame( M_d)
  names(M_d) <- c(paste("x.lag", kMax:0, sep = ""), 
                  paste("y.lag", kMax:0, sep = ""))
  
  # regression : x_0 ~ x_1 +...+ x_jj + y_0 +...+ y_jj, and selection of the maximun order (p1) 
  # using BIC criteria 
  # the jj value is between 0 and k value. (see inputs in the header)
  
  bic1 <- data.frame()
  
  for (jj in kinf:(kMax)){
    b1 <- BIC(model <- lm(x.lag0 ~ . , data = M_d[, c((kMax+1-jj):(kMax+1), 
                                                      (2*(kMax+1)-jj):(2*(kMax+1)))]))
    bic1 <- rbind(bic1, c(jj, b1))
  }
  
  names(bic1) <- c("k", "BIC")
  p1 <- bic1[which.min(bic1$BIC), "k"]
  
  # if in the linear regressions the maximun order p1 is equal to k we stop the function.  
  # otherwise (p1 < k) the next part determines the maximum order using the regressions
  # y_0 ~ y_1 +...+ y_jj + x_0 +...+ x_jj, 
  # and selection of the order p2 using BIC criteria
  
  
  if(p1 < kMax){
    bic2 <- data.frame()
    
    for (jj in kinf:(kMax)){
      b2 <- BIC(model <- lm(y.lag0 ~ . , data = M_d[, c((kMax+1-jj):(kMax+1), 
                                                        (2*(kMax+1)-jj):(2*(kMax+1)))]))
      bic2 <- rbind(bic2, c(jj, b2))
    }
    
    names(bic2) <- c("k", "BIC")
    p2 <- bic1[which.min(bic2$BIC), "k"]
    p  <- max(p1, p2)
    
  }
  
  if(!exists("p")) p <- p1
  
  return(p)
}


##############################################################################################################################################

#.............................................................................
#
# WithinDispersion.R
# 
# Date : 19-04-2018
# 
# Author:              Carolina Gamboa. 
# 
# Description:   Within dispersion measures. Expression (2) at Tibshirani et al (2001).
# 
# Input :                    DistanceMatrix : square matrix of GCC distances
#                            Clusters : Matriz de asignación de grupos 
#                            nClus: number of groups 
# 
# Output:                   Within dispersion measure 
# 
#.............................................................................

WithinDispersion = function(DistanceMatrix, Clusters, nClus){
  RW <- NULL
  for (k in 1:nClus){
    
    D <- NULL
    n <- NULL
    for (r in 1:k){
      Indexes <- which(Clusters[,k] == r)
      D[r] <- sum(sum(DistanceMatrix[Indexes,Indexes]))
      n[r] = 2*sum(length(Indexes));
    }
    RW[k] <- sum(D/n)
  }
  return(RW)
}

###########################################################################################################

#.............................................................................
#
# GAPdistance.R
# 
# Date : 19-04-2018
# 
# Author:              Carolina Gamboa. 
# 
# Description:  this function computes the gap and the number of groups using the gap statistics
# 
# Input :                    DistanceMatrix : square matrix of GCC distances
#                            Clusters : Matriz de asignación de grupos 
#                            B: number of iterations for the bootstrap  
# 
# Output:                   nClus number of groups 
# 
#.............................................................................

GAPdistance <- function(DistanceMatrix, Clusters, B){
  
  # multidimentional scaling
  
  N <- dim(Clusters)[1]
  nClus <- dim(Clusters)[2]
  
  # % Within dispersion measures at the observed data.
  W <-  WithinDispersion(DistanceMatrix, Clusters, nClus);
  
  mds <- cmdscale(DistanceMatrix,eig=TRUE)
  eps <- 2^(-52)
  f <-  sum(mds$eig > eps^(1/4))
  mds <- cmdscale(DistanceMatrix,eig=TRUE, k = f)
  
  Xmatrix <- mds$points
  
  # % PCA reference feature space.
  svd <- svd(Xmatrix); U <- svd$u; V <- svd$v; D <- svd$d;
  Zmatrix = Xmatrix %*% V;
  Zmin = apply(Zmatrix, 2, min);
  Zmax = apply(Zmatrix, 2, max);
  
  
  #% Within dispersion measures at the reference feature space.
  Wstar = matrix(ncol = B, nrow = nClus)
  for (b in 1:B){
  
    for (ff in 1:f){
      Zmatrix[ ,ff] = runif(N, Zmin[ff], Zmax[ff]);
    }
    
      Zmatrix <- Zmatrix %*% t(V);
      ZDistanceMatrix = (dist(Zmatrix));
      L = hclust(dist(Zmatrix), method = "single");
      ZClusters = cutree(L, k = 1:nClus);
      ZDistanceMatrix <- as.matrix(ZDistanceMatrix)
      Wstar[,b] = WithinDispersion(ZDistanceMatrix, ZClusters, nClus);
  
  }
  
  
  logWmean = apply(log(Wstar), 1, mean);
  logWstd  = apply(log(Wstar), 1, sd)*sqrt(1 + 1/B);
  
  
  GAPstat = logWmean - log(W);
  
  WhoseK = GAPstat[1:nClus-1] - GAPstat[2:nClus] + logWstd[2:nClus];
  
  R <-  min(which(WhoseK >= 0))
  return(R)

}
  

###########################################################################################################

#.............................................................................
#
# graphMatrix.R
# 
# Date : 19-04-2018
# 
# Author:              Carolina Gamboa. 
# 
# Description:  this function presents a graph with the distribution of gcc distances 
#               for each couple of groups
# 
# Input :                    DistanceMatrix : square matrix of GCC distances
#                            Clusters : Matriz de asignación de grupos 
# 
# Output:                   graph square
# 
#.............................................................................


graphMatrix <- function(distMatrix, clusters){
  
  colnames(distMatrix) <- 1:nrow(distMatrix)
  rownames(distMatrix) <- 1:nrow(distMatrix)
  
  nClus <- length(unique(clusters))
  par(mfrow = c(nClus, nClus))
  for(i in c(1:max(clusters), 0)){
    for(j in c(1:max(clusters), 0)){
      subset <- which(clusters%in%c(i, j))
      sampleDM <- as.dist(distMatrix[subset, subset])
      
      if(i != j){ 
        mainPlot <- paste("Density of distances for clusters", i, "and",j)
      } else {
        mainPlot <- paste("Density of distances for cluster", j)
      }
      plot(density(sampleDM), main = mainPlot, 
           xlab = paste("median dist ", round(median(sampleDM), 4), sep =""))
      }
  }
  
  
}

#######################################################################################################################################################


############################################################################################################################################################
#.............................................................................
#
# silhouetteNClus.R
# 
# Date  : 19-04-2018
# 
# Author:              Carolina Gamboa. 
# 
# Description:  Find the number of cluster by the Silhouette statistics

# Input :                 nClus : Maximum number of groups
#                            dis : Matrix of GCC distances  
#                            method: Hierarchical method "single", "average"  
# 
# Output:                  nClus : Number of groups
#                         list : Silhouette statistics for each value of nclus 
# 
#.............................................................................

#library

library(cluster)

silhouetteNClus <- function(nClus, distanceMatrix, method){
  
  if(class(distanceMatrix)!="dist") distanceMatrix <- as.dist(distanceMatrix)
  
  silIndex <- data.frame(nClus = 1, silIndex = 0)
  
  if(method == "complete"){
    
    
    hclust.dist <- hclust(distanceMatrix, method="complete")
    Cl.hclust <- cutree(hclust.dist, 2:nClus)
    
    for (jj in 2:nClus) {
      coef <- silhouette(Cl.hclust[, jj-1], distanceMatrix)
      jjSilIndex <- mean(coef[, "sil_width"])
      silIndex <- rbind(silIndex, data.frame(nClus = jj, silIndex = jjSilIndex))
      
    }
    
    
  }
  
  if(method == "average"){
    
    
    hclust.dist <- hclust(distanceMatrix, method="average")
    Cl.hclust <- cutree(hclust.dist, 2:nClus)
    
    for (jj in 2:nClus) {
      coef <- silhouette(Cl.hclust[, jj-1], distanceMatrix)
      jjSilIndex <- mean(coef[, "sil_width"])
      silIndex <- rbind(silIndex, data.frame(nClus = jj, silIndex = jjSilIndex))
      
    }
    
    
  }
  
  
  
  if(method == "single"){
    
    
    hclust.dist <- hclust(distanceMatrix, method="single")
    Cl.hclust <- cutree(hclust.dist, 2:nClus)
    
    for (jj in 2:nClus) {
      coef <- silhouette(Cl.hclust[, jj-1], distanceMatrix)
      jjSilIndex <- mean(coef[, "sil_width"])
      silIndex <- rbind(silIndex, data.frame(nClus = jj, silIndex = jjSilIndex))
      
    }
    
    
  }
  maxPos <- which(silIndex[, "silIndex"]==max(silIndex[, "silIndex"]))
  cluster <- list(nClus = silIndex[maxPos, "nClus"], coef = silIndex)
  return(cluster$nClus)
  
}

    
  
  
    
######################################################################################################################################################
#########################################################################################################################################################