#' @title Sparse PCA 
#' @export
#' 
#' @description This function performs sparse principal component analysis (PCA)  
#' using a block optimisation algorithm or an iterative deflation algorithm.
#' 
#' @param A a numerical data matrix of size n by p (observations by variables)
#' @param m the number of components
#' @param lambda a numerical vector of size m providing the reduced sparsity 
#' parameters (in relative value with respect to the theoretical upper bound). 
#' Each reduced sparsity parameter is a value between 0 and 1.
#' @param block either 0 or 1. block==0 means that deflation is used if more 
#' than one component. A block optimisation algorithm is otherwise used
#' that computes m components at once. By default, block=1.
#' @param mu numerical vector of size m with the mu parameters (required for the block algorithms only).
#' By default, mu_j=1/j.
#' @param center a logical value indicating whether the variables should be 
#' shifted to be zero centered.
#' @param scale a logical value indicating whether the variables should be 
#' scaled to have unit variance.
#' @param iter_max maximum number of admissible iterations.
#' @param epsilon accuracy of the stopping criterion. 
#' 
#' @return \item{Z}{a p by m numerical matrix with the m sparse loading vectors} 
#' @return \item{Y}{a n by m numerical matrix with the m principal components} 
#' @return \item{B}{the numerical data matrix centered (if center=TRUE) and scaled 
#' (if scale=TRUE)}
#' 
#' @details This function implements an optimal projected variance sparse
#' block PCA algorithm and a deflation algorithm applying the block algorithm with 
#' one single component iteratively to each deflated data matrix. 
#' 
#' The block algorithm uses a numerical vector of parameters \code{mu} usually 
#' chosen either striclty decreasing (mu_j=1/j) or all equal (mu_j=1 for all j).
#' Striclty decreasing parameters relieve the underdetermination which happens 
#' in some situations and drives to a solution close to the PCA solution.
#' 
#' The principal components are defined by Y=BZ  where B is the centered (if 
#' center=TRUE) and scaled (if scale=TRUE) data matrix and where Z is the sparse loading matrix. 
#' 
#'@references 
#'M. Chavent and G. Chavent, Optimal projected variance group-sparse block PCA, 
#'submitted, 2020.
#'@references
#'M. Journee, Y. Nesterov, P. Richtarik, and R. Sepulchre. Generalized power 
#'method for sparse principal component analysis. Journal of Machine Learning
#'Research, 11:517-553, 2010.
#'
#'@seealso \code{\link{groupsparsePCA}}, \code{\link{pev}}, 
#'\code{\link{explainedVar}}
#'
#' @examples 
#' # Simulated data
#'  v1 <- c(1,1,1,1,0,0,0,0,0.9,0.9)
#'  v2 <- c(0,0,0,0,1,1,1,1,-0.3,0.3)
#'  valp <- c(200,100,50,50,6,5,4,3,2,1)
#'  A <- simuPCA(50,cbind(v1,v2),valp,seed=1)
#'  
#'  # Three sparse PCA algorithms
#'  Z <- sparsePCA(A,2,c(0.5,0.5),block=0)$Z #deflation
#'  Z <- sparsePCA(A,2,c(0.5,0.5),block=1)$Z #block different mu
#'  Z <- sparsePCA(A,2,c(0.5,0.5),block=1,mu=c(1,1))$Z #block same mu
#'  
#'  # Example of the protein data
#'  data("protein")
#'  Z <- sparsePCA(protein,2,c(0.5,0.5))$Z #block different mu
#'  
#'
sparsePCA <- function(A,m,lambda, block=1, mu=1/1:m,
                      center=TRUE,scale=TRUE,iter_max=1000,epsilon=0.0001  )
{
  n <- nrow(A)
  p <- ncol(A)
  if (m > min(n-1,p))
     stop("m must be smaller than rank(A)",call. = FALSE)
  if (length(lambda)!=m) 
    stop("lambda must be a vector of size m",call. = FALSE)
  if ((max(lambda) >1) || (max(lambda) < 0))
    stop("Values in lambda must be in [0,1]",call. = FALSE)
  
  norm2 <- function(x) sqrt(sum(x^2))
  polar <- function(x)
  {
    obj <- svd(x)
    obj$u %*% t(obj$v)
  }
  
  Z <- matrix(0,p,m)
  rownames(Z) <- colnames(A)
  A <- as.matrix(A)

  
  if (center==TRUE)
    A <- scale(A,center=TRUE, scale=FALSE) #center data
  if (scale==TRUE)
    A <- scale(A,center=FALSE,scale=TRUE)*sqrt(n/(n-1)) #scale data
 
  if ((m==1) || (m>1 && block==0)) # deflation is used if m>1
  {
    B <- A
    gamma <- rep(NA,m)
    X <- matrix(NA,n,m)
    for (comp in 1:m)   #loop on the components
    {
      i_max <- which.max(apply(B,2,norm2)) 
      gamma_max <- norm2(B[,i_max])
      gamma[comp] <- lambda[comp]*gamma_max
      x <-svd(B)$u[,1]     #initialization
      f <- matrix(0,iter_max,1)
      iter <- 1
      repeat
      {
        Bx <- t(B)%*%x
        tresh <- sign(Bx)*(abs(Bx)>=gamma[comp])*(abs(Bx)-gamma[comp])
        f[iter] <- sum(tresh^2)  #cost function
        if (f[iter]==0) #the sparsity parameter is too high: all entries of the loading vector are zero
          break
        grad <- B %*% tresh
        x <- grad/norm2(grad)
        if (((iter>=2) && ((f[iter]-f[iter-1])/f[iter-1] < epsilon)) || (iter >= iter_max)) #stopping criterion
        {
          if (iter >= iter_max) print("Maximum number of iterations reached")
          break
        }
        iter <- iter+1
      }
      Bx <- t(B)%*%x
      #pattern <- abs(Bx)-gamma >0
      z <- sign(Bx)*(abs(Bx)>=gamma[comp])*(abs(Bx)-gamma[comp]) 
      if (max(abs(z))>0)
        z <- z/norm2(z)
      #if (post==TRUE) z <- pattern_filling(B,z)
      y <- B%*%z  
      B <- B-y%*%t(z)
      #B <- scale(B,center=TRUE,scale=FALSE) # center the matrix of residuals
      #B <- scale(B)*sqrt(n/(n-1)) # standardize the matrix of residuals
      Z[,comp]<-z
      X[,comp] <- y/norm2(y)
    }
  }
  if (m>1 && block==1) # block algorithm
  {
    e <- svd(A)
    d <- e$d[1:m]
    u <- e$u[,1:m]
    i_max <- which.max(apply(A,2,norm2)) 
    gamma_max <- norm2(A[,i_max])
    gamma <- lambda*gamma_max*d/d[1]
    x <-u     #initialization 
    f <- matrix(0,iter_max,1)
    iter <- 1
    repeat
    {
      Ax <- t(A)%*%x
      tresh <- sign(Ax)*sweep(abs(Ax),2,gamma,">=")*sweep(abs(Ax),2,gamma,"-")
      #f[iter] <-Tr(crossprod(tresh)%*%diag(mu^2)) #cost function
      f[iter] <- sum(apply(tresh,2,norm2)^2*mu^2) #cost function
      if (f[iter]==0) #the sparsity parameter is too high: all entries of the loading matrix are zero
        break
      grad <- 2*A %*% tresh%*%diag(mu^2)
      x <- polar(grad)
      if (((iter>=2) && ((f[iter]-f[iter-1])/f[iter-1] < epsilon)) || (iter >= iter_max)) #stopping criterion
      {
        if (iter >= iter_max) print("Maximum number of iterations reached")
        break
      }
      iter <- iter+1
    }
    Ax <- t(A)%*%x
    z <- sign(Ax)*sweep(abs(Ax),2,gamma,">=")*sweep(abs(Ax),2,gamma,"-")
    for (j in 1:m)
    {
      if (max(abs(z[,j]))>0)    
        Z[,j] <- z[,j]/norm2(z[,j])
    }
    #if (post==TRUE) Z <- pattern_filling(A,Z,mu)
    X=x
  }
  Y <- A %*% Z
  
  res <- list(Z=Z,Y=Y,gamma=gamma,B=A,X=X)
  class(res) <- "sparsePCA"
  return(res)
}


