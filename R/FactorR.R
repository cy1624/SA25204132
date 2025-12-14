Sigma_k <- function(data, k){
  ybar <- colMeans(data)
  n <- dim(data)[1]
  p <- dim(data)[2]
  Sigma <- cov(data[1:(n-k),],data[(k+1):n,])
  return(Sigma)
}



#' @title Calculating the numbers of factors
#' @description This function calculates the numbers of both common factors and specific factors. See more in \url{https://doi.org/10.1080/01621459.2023.2183132}.
#' @param data a data frame or matrix
#' @param k0 a pre-specified integer used in calculating covariance matrix. The exact value won't affect the result much and it could be proper to choose it as 5(the default value) in practice.
#' @param j0 a pre-specified integer representing the maximum of the number of factors. It must be less than the number of variables. The default value is one forth of the number of variables(rounded down to the nearest integer).
#' @return a list consists of \code{r0}, \code{r}, \code{R} and \code{tao}. \code{r0} is the number of common factors. \code{r} is the number of specific factors. \code{R} is all the ratios of eigenvalue sums. \code{tao} is the indices of local maximums of R, which is sorted into descending order by the corresponding value of R.
#' @examples
#' \dontrun{
#' data(CSI300_2)
#' nums <- FactorNum(CSI300_2)   
#' cat(c(nums$r0, nums$r))
#' }
#' @export
FactorNum <- function(data, k0=5, j0=ceiling(dim(data)[2]/4)){
  n <- dim(data)[1]
  p <- dim(data)[2]
  if(is.data.frame(data)) data <- as.matrix(data)
  EigenValues <- sapply(0:(k0+1), function(k){
    eigen(Sigma_k(data,k)%*%t(Sigma_k(data,k)))$values[1:(j0+1)]
  })
  EigenSums <- rowSums(EigenValues)
  R <- EigenSums[1:j0]/EigenSums[2:(j0+1)]
  R <- c(1,R)
  LocalMax_index <- c()
  for(i in 2:j0){
    if(R[i]>R[i-1]&R[i]>R[i+1]) LocalMax_index <- c(LocalMax_index,i)
  }
  tao <- LocalMax_index[order(-R[LocalMax_index])] - 1
  tao1 <- tao[1:2]
  return(list(r0=min(tao1), r=max(tao1)-min(tao1), R=R[-1], tao=tao))
}

#' @title Calculating the loading matrices of factors
#' @description This function calculates the loading matrices of both common factors and specific factors. See more in \url{https://doi.org/10.1080/01621459.2023.2183132}.
#' @param data a data frame or matrix
#' @param k0 a pre-specified integer used in calculating the numbers of factors. It works only when r is \code{NULL}.
#' @param j0 a pre-specified integer used in calculating the numbers of factors. It works only when r is \code{NULL}.
#' @param r the numbers of factors. The first value is the numbers of common factors and the second value is the numbers of specific factors. When r is \code{NULL}, the function \code{FactorNum} will be used to calculate the numbers of the factors.
#' @return a list consists of \code{Loadings_CommonFactor} and \code{Loadings_SpecificFactor}, which represent the loadings of common factors and specific factors, respectively.
#' @examples
#' \dontrun{
#' data(CSI300_2)
#' Loadings_Factor(CSI300_2, r=c(1,13))
#' }
#' @export
Loadings_Factor <- function(data, k0=5, j0=ceiling(dim(data)[2]/4), r=NULL){
  n <- dim(data)[1]
  p <- dim(data)[2]
  if(is.data.frame(data)) data <- as.matrix(data)
  M1 <- matrix(0, p, p)
  for(i in 0:k0){
    M1 <- M1 + Sigma_k(data, i)%*%t(Sigma_k(data, i))
  }
  if(is.null(r)){
    r <- FactorNum(data,k0,j0)
    r0 <- r$r0
    r1 <- r$r
  }else{
    r0 <- r[1]
    r1 <- r[2]
  }
  
  A <- eigen(M1)$vectors[,1:r0]
  data1 <- data %*% (diag(1,p,p)-A%*%t(A))
  M2 <- matrix(0, p, p)
  for(i in 0:k0){
    M2 <- M2 + Sigma_k(data1, i)%*%t(Sigma_k(data1, i))
  }
  B <- eigen(M2)$vectors[,1:r1]
  return(list(Loadings_CommonFactor=A, 
              Loadings_SpecificFactor=B))
}

#' @title Preparation for clustering
#' @description This function calculates the the components not belonging to any clusters, the number of clusters and the correlation matrix based on the loading of specific factors, which can be used in k-means clustering. See more in \url{https://doi.org/10.1080/01621459.2023.2183132}.
#' @param data a data frame or matrix
#' @param loadings a list containing \code{Loadings_SpecificFactor}. It should be like \code{list(Loadings_SpecificFactor=...)}. If \code{NULL}, the function \code{Loadings_Factor} will be used to calculate the loadings of the factors.
#' @param k0 a pre-specified integer used in calculating the numbers of factors. It works only when loadings and r is \code{NULL}.
#' @param j0 a pre-specified integer used in calculating the numbers of factors. It works only when loadings and r is \code{NULL}.
#' @param r the numbers of factors. If \code{NULL}, the function \code{FactorNum} will be used to calculate the numbers of the factors.
#' @param omega the threshold deciding the components not belonging to any clusters. If \code{NULL}, it will be set to a proper number according to the data and loadings.
#' @return a list consists of \code{NoCluster_Index}, \code{d} and \code{R}. \code{NoCluster_Index} is the indices of the components not belonging to any clusters. \code{d} is the number of clusters. \code{R} is the correlation matrix based on the loading of specific factors, which can be used in k-means clustering.
#' @examples
#' \dontrun{
#' data(CSI300_2)
#' Clusters(CSI300_2)$d
#' }
#' @export
Clusters <- function(data, loadings=NULL,  k0=5, j0=ceiling(dim(data)[2]/4), r=NULL, omega=NULL){
  n <- dim(data)[1]
  p <- dim(data)[2]
  if(is.data.frame(data)) data <- as.matrix(data)
  if(is.null(loadings)){
    loadings <- Loadings_Factor(data,k0,j0,r)
  }
  
  B <- loadings$Loadings_SpecificFactor
  r <- dim(B)[2]
  
  if(is.null(omega)){
    omega <- sqrt(r/(p*log(p)))
  }
  
  norms <- apply(B, 1, FUN = function(x) norm(x,type = "2"))
  NoCluster_Index <- which(apply(B, 1, FUN = function(x) norm(x,type = "2"))<omega)
  d <- sum(eigen(abs(B%*%t(B)))$values > 1-1/log(n))
  p0 <- p-length(NoCluster_Index)
  rownames(B) <- colnames(data)
  F <- B[-NoCluster_Index,]
  R <- matrix(0, p0, p0)
  for(i in 1:p0){
    for(j in 1:p0){
      R[i,j] <- abs(t(F[i,])%*%F[j,])/sqrt((t(F[i,])%*%F[i,])*(t(F[j,])%*%F[j,]))
    }
  }
  
  return(list(NoCluster_Index=NoCluster_Index,d=d,R=R))
}





#' @import Rcpp 
#' @import stats 
#' @import RcppArmadillo
#' @import rbenchmark
#' @import magrittr
#' @import snowfall
#' @import boot
#' @import lpSolve
#' @import microbenchmark
#' @import DAAG 
#' @import bootstrap
#' @importFrom gtools permutations
#' @useDynLib SA25204132
NULL

#' @title Daily closing prices of CSI300 components
#' @name CSI300_1
#' @description This dataset comprises the daily closing prices of selected constituent stocks of the CSI 300 Index from January 1, 2018 to December 31, 2023, with missing values subjected to interpolation. The column names correspond to the stock codes. 
#' @examples
#' \dontrun{
#' data(CSI300_1)
#' head(CSI300_1[,1:15])
#' }
NULL

#' @title Logarithmic return of CSI300 components
#' @name CSI300_2
#' @description This dataset comprises the logarithmic return of the daily closing prices of selected constituent stocks of the CSI 300 Index from January 1, 2018 to December 31, 2023, which is derived from CSI300_1.
#' @examples
#' \dontrun{
#' data(CSI300_2)
#' head(CSI300_2[,1:15])
#' }
NULL