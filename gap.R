

#.............................................................................
#
# gap.R
# 
# Fecha de Creación : 19-04-2018
# 
# Autor:              Carolina Gamboa. 
# 
# Descripción:  Esta función calcula el indice gap y determina la cantidad 
#               de grupos optima
# 
# Entradas: Tres parámetros. DistanceMatrix : matriz de distancia
#                            Clusters : Matriz de asignación de grupos 
#                            B: cantidad de iteraciones bootstrap  
# 
# Output: Dos parámetros. K : Cantidad de grupos
# 
#.............................................................................



WithinDispersion = function(DistanceMatrix, Clusters, K){
  RW <- NULL
  for (k in 1:K){
    
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




GAPdistance <- function(DistanceMatrix, Clusters, B){
  
  # escalamiento multidimensional
  
  N <- dim(Clusters)[1]
  K <- dim(Clusters)[2]
  
  # % Within dispersion measures at the observed data.
  W <-  WithinDispersion(DistanceMatrix, Clusters, K);
  
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
  Wstar = matrix(ncol = B, nrow = K)
  for (b in 1:B){
  
    for (ff in 1:f){
      Zmatrix[ ,ff] = runif(N, Zmin[ff], Zmax[ff]);
    }
    
      Zmatrix <- Zmatrix %*% t(V);
      ZDistanceMatrix = (dist(Zmatrix));
      L = hclust(dist(Zmatrix), method = "single");
      ZClusters = cutree(L, k = 1:K);
      ZDistanceMatrix <- as.matrix(ZDistanceMatrix)
      Wstar[,b] = WithinDispersion(ZDistanceMatrix, ZClusters, K);
  
  }
  
  
  logWmean = apply(log(Wstar), 1, mean);
  logWstd  = apply(log(Wstar), 1, sd)*sqrt(1 + 1/B);
  
  
  GAPstat = logWmean - log(W);
  
  WhoseK = GAPstat[1:K-1] - GAPstat[2:K] + logWstd[2:K];
  
  R <-  min(which(WhoseK >= 0))
  return(R)

}
  



# % Within dispersion measures. Expression (2) at Tibshirani et al (2001).

# 
# X <- scale(iris[, -5])
# N <- nrow(X)
# euc.distance <- dist(X)
# linka <- hclust(euc.distance, method = "single")
# plot(linka, labels = FALSE, hang = -1)
# maxcluster <- cutree(linka, k = 1:10)
# 
# WithinDispersion(as.matrix(euc.distance), cluster, 10)
# GAPdistance(as.matrix(euc.distance), cluster, 100)



