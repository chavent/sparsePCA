#' @title Sparse PCA 
#' @description Sparse principal component analysis (PCA). Three algorithms are proposed for the determination of the
#' sparse loadings : on based on a deflation procedure and two based on block optimisation.
#' @param A a n times p  numerical data matrix
#' @param m number of components
#' @param lambda a vector of dimension m with reduced sparsity parameters (in relative value with respect to the theoretical upper bound). 
#' Each reduced sparsity parameter is a value between 0 and 1.
#' @param block either 0 or 1. block==0 means that deflation is used if more than one component need to be computed.
#' A block algorithm is otherwise used, that computes m components at once. By default, block=1.
#' @param mu vector of dimension m with the mu parameters (required for the block algorithms only).
#' By default, mu_j=1/j.
#' @param center a logical value indicating whether the variables should be shifted to be zero centered.
#' @param scale a logical value indicating whether the variables should be scaled to have unit variance.
#' @param iter_max maximum number of admissible iterations.
#' @param epsilon accuracy of the stopping criterion. 
#' @param post If TRUE the loadings are obtained with post-processing.
#' @return \item{Z}{the p times m matrix that contains the m sparse loading vectors} 
#' @return \item{Y}{the n times m matrix that contains the m principal components} 
#' @return \item{B}{the data matrix centered (if center=TRUE) and scaled (if scale=TRUE)}
#' @details The principal components are given by Y=BZ  where B is the matrix A which has been centered (if scale=FALSE) or standardized (if scale=TRUE) and Z is the matrix of the sparse loading vectors. 
#'@references 
#'\itemize{
#'\item M. Chavent and G. Chavent, Group-sparse block PCA and explained variance, arXiv:1705.00461
#'\item M. Journee and al., Generalized Power Method for Sparse Principal Component Analysis, Journal of Machine Learning Research 11 (2010) 517-553.
#'\item H. Shen and J.Z. Huang, Journal of Multivariate Analysis 99 (2008) 1015-1034.
#'}
#' @export
#' @examples 
#' # Example from Shen & Huang 2008
#'  v1 <- c(1,1,1,1,0,0,0,0,0.9,0.9)
#'  v2 <- c(0,0,0,0,1,1,1,1,-0.3,0.3)
#'  valp <- c(200,100,50,50,6,5,4,3,2,1)
#'  n <- 50
#'  A <- simuPCA(n,cbind(v1,v2),valp,seed=1)
#'  # Three algorithms
#'  Z <- sparsePCA(A,2,c(0.5,0.5),block=0)$Z #deflation
#'  Z <- sparsePCA(A,2,c(0.5,0.5),block=1)$Z #block different mu
#'  Z <- sparsePCA(A,2,c(0.5,0.5),block=1,mu=c(1,1))$Z #block same mu
#'  
#'  # Example of the protein data
#'  data("protein")
#'  Z <- sparsePCA(protein,2,c(0.5,0.5))$Z #block different mu
#'  
#'@seealso \code{\link{groupsparsePCA}}, \code{\link{pev}}
#'
sparsePCA <- function(A,m,lambda,block=1, mu=1/1:m,post=FALSE,center=TRUE,scale=TRUE,iter_max=1000,epsilon=0.0001  )
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
  
  iter_max <- 1000   # maximum number of admissible iterations
  epsilon <- 0.0001  # accuracy of the stopping criterion 

  Z <- matrix(0,p,m)
  rownames(Z) <- colnames(A)
  A <- as.matrix(A)
  gamma <- rep(NA,m)
  
  if (center==TRUE)
    A <- scale(A,center=TRUE, scale=FALSE) #center data
  if (scale==TRUE)
    A <- scale(A,center=FALSE,scale=TRUE)*sqrt(n/(n-1)) #scale data
 
  if ((m==1) || (m>1 && block==0)) # deflation is used if m>1
  {
    B <- A
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
      #B <- scale(B<)*sqrt(n/(n-1)) # standardize the matrix of residuals
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


#' Ztrue data
#' @format A data  matrix with 20 rows and 4 columns. 
#' @source Chavent & Chavent (2017).
#' @description The four columns are specified leading eigenvectors of a covariance matrix used to
#' simulate data from a group-sparse PCA model : the 20 variables are organized in
#' 5 groups of 4 variables. 
#' @name Ztrue
NULL

#' Protein data
#' @format A data matrix with 25 rows (the European countries) and 9 columns (the food groups) 
#' @source Originated by A. Weber and cited in Hand et al., A Handbook of Small Data Sets, (1994, p. 297).
#' @description The data measure the amount of protein consumed for nine food groups in 
#' 25 European countries. The nine food groups are red meat (RedMeat), white meat (WhiteMeat), 
#' eggs (Eggs), milk (Milk), fish (Fish), cereal (Cereal), starch (Starch), nuts (Nuts), and fruits and vegetables (FruitVeg).
#' @name protein 
NULL