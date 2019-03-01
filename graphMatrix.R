

graphMatrix <- function(distMatrix, clusters){
  
  colnames(distMatrix) <- 1:nrow(distMatrix)
  rownames(distMatrix) <- 1:nrow(distMatrix)
  
  nClusters <- length(unique(clusters))
  par(mfrow = c(nClusters, nClusters))
  for(i in c(1:max(clusters), 0)){
    for(j in c(1:max(clusters), 0)){
      subset <- which(clusters%in%c(i, j))
      sampleDM <- as.dist(distMatrix[subset, subset])
      
      if(i != j){ 
        mainPlot <- paste("Density of distances for clusters", i, "and",j)
      } else {
        mainPlot <- paste("Density of distances for cluster", j)
      }
      plot(density(sampleDM), main = mainPlot)
      }
  }
  
  

  
}