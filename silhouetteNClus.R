#.............................................................................
#
# silhouetteNClus.R
# 
# Fecha de Creación : 19-04-2018
# 
# Autor:              Carolina Gamboa. 
# 
# Descripción:  Esta función calcula el indice de silueta para determinadas 
#               cantidades de grupos, y determina la cantidad de grupos optima
#               usando el estadístico de silueta. 
# 
# Entradas: Tres parámetros. K : Cantidad máxima de grupos a analizar
#                            dis : Matriz de distancias. 
#                            method: Método Jerárquico "single", "average"  
# 
# Output: Dos parámetros. K : Cantidad de grupos
#                         list : valor del índice de silueta para cada k 
# 
#.............................................................................

#library

library(cluster)

silhouetteNClus <- function(K, dis, method){
  
  if(class(dis)!="dist") dis <- as.dist(dis)
  
  silIndex <- data.frame(k = 1, silIndex = 0)
  
  if(method == "complete"){
 
      
      hclust.dist <- hclust(dis, method="complete")
      Cl.hclust <- cutree(hclust.dist, 2:K)
      
      for (jj in 2:(K)) {
      coef <- silhouette(Cl.hclust[, jj-1], dis)
      jjSilIndex <- mean(coef[, "sil_width"])
      silIndex <- rbind(silIndex, data.frame(k = jj, silIndex = jjSilIndex))
      
    }
    
    
  }

  if(method == "average"){
    
    
    hclust.dist <- hclust(dis, method="average")
    Cl.hclust <- cutree(hclust.dist, 2:K)
    
    for (jj in 2:(K)) {
      coef <- silhouette(Cl.hclust[, jj-1], dis)
      jjSilIndex <- mean(coef[, "sil_width"])
      silIndex <- rbind(silIndex, data.frame(k = jj, silIndex = jjSilIndex))
      
    }
    
    
  }
  
  
  
  if(method == "single"){
    
    
    hclust.dist <- hclust(dis, method="single")
    Cl.hclust <- cutree(hclust.dist, 2:K)
    
    for (jj in 2:(K)) {
      coef <- silhouette(Cl.hclust[, jj-1], dis)
      jjSilIndex <- mean(coef[, "sil_width"])
      silIndex <- rbind(silIndex, data.frame(k = jj, silIndex = jjSilIndex))
      
    }
    
    
  }
  maxPos <- which(silIndex[, "silIndex"]==max(silIndex[, "silIndex"]))
  cluster <- list(K = silIndex[maxPos, "k"], coef = silIndex)
  return(cluster$K)

  }
  