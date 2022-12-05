#'@title  Simulation from a PCA model
#'@export
#'
#'@description This function  simulates data from a PCA model. 
#'The data matrix is drawn from a centered multivariate normal 
#'distribution where the covariance matrix is constructed 
#'knowing its p eigenvalues and its m < p first eigenvectors 
#'(see Chavent & Chavent 2020).
#'
#'@param n the number of observations 
#'@param Ztrue a numerical matrix of size p by m with the underlying loading 
#'vectors (m first eigenvectors of the underlying covariance
#'matrix)
#'@param valp a numerical vector of size p with the eigenvalues of 
#'the underlying covariance matrix
#'@param seed a seed value to replicate the simulation. By default, seed=FALSE
#'
#'@return Returns a numerical data matrix of size n by p
#'
#'@references 
#'M. Chavent and G. Chavent, Optimal projected variance group-sparse block PCA, 
#'submitted, 2020.
#'@references
#'M. Journee, Y. Nesterov, P. Richtarik, and R. Sepulchre. Generalized power 
#'method for sparse principal component analysis. Journal of Machine Learning
#'Research, 11:517-553, 2010.
#'@references H. Shen and J.Z. Huang. Sparse principal component analysis via 
#' regularized low rank matrix approximation. Journal of multivariate analysis, 
#' 99(6):1015-1034, 2008.
#'
#'@seealso \code{\link{sparsePCA}}, \code{\link{groupsparsePCA}}
#'
#'@examples
#' # Example from Journee & al. 2010
#'  v1 <- c(rep(1/sqrt(10),10),rep(0,490))
#'  v2 <- c(rep(0,10), rep(1/sqrt(10),10),rep(0,480))
#'  valp <- c(400,300,rep(1,498))
#'  n <- 100
#'  A <- simuPCA(n,cbind(v1,v2),valp,seed=1)
#'  svd(A)$d[1:3]^2 #eigenvalues of the empirical covariance matrix.
#'  
#' # Example from Shen & Huang 2008
#'  v1 <- c(1,1,1,1,0,0,0,0,0.9,0.9)
#'  v2 <- c(0,0,0,0,1,1,1,1,-0.3,0.3)
#'  valp <- c(200,100,50,50,6,5,4,3,2,1)
#'  n <- 100
#'  A <- simuPCA(n,cbind(v1,v2),valp,seed=1)
#'  svd(A)$d^2 #eigenvalues of the empirical covariance matrix (times n).
#'  
#' # Example from Chavent & Chavent 2020
#'  valp <- c(c(200,180,150,130),rep(1,16))
#'  data("Ztrue")
#'  n <- 100
#'  A <-simuPCA(n,Ztrue,valp,seed=1)
#'  svd(A)$d^2 #eigenvalues of the empirical covariance matrix.
#'  

simuPCA <- function(n=30, Ztrue, valp, seed=FALSE)
{
  norm2 <- function(x) sqrt(sum(x^2))
  p <- nrow(Ztrue)
  m <- ncol(Ztrue)
  if (m > p) stop("the number of loadings (number of columns in Ztrue) 
                  must be smaller of equal to the number of variables (rowns)",call.=FALSE)
  for (j in 1:ncol(Ztrue))
        Ztrue[,j] <- Ztrue[,j]/norm2(Ztrue[,j])
  if (seed != FALSE) set.seed(seed)

  if (length(valp)!=p)
    stop("the length of valp must be equal to the number of rows in Ztrue",call. = FALSE)
  U <- matrix(stats::runif(p*(p-m)),p,(p-m))
  V <- qr.Q(qr(cbind(Ztrue,U)))
  C <- V%*%diag(valp)%*%t(V)
  if (seed != FALSE) set.seed(seed)
  A <- mvtnorm::rmvnorm(n,sigma=C) /sqrt(n)
  return(A)
}

