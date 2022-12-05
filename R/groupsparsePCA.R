#' @title Group-sparse PCA 
#' @export
#' 
#' @description This function implements group-sparse principal component 
#' analysis (PCA) using a block optimisation algorithm or an iterative 
#' deflation algorithm. It generalizes the block optimisation approach of Journee 
#' et al. (2010) to group sparsity.
#' 
#' @param A a numerical data matrix of size n by p (observations by variables)
#' @param m number of components
#' @param lambda  a numerical vector of size m providing reduced sparsity parameters 
#' (in relative value with respect to the theoretical upper bound). 
#' Each reduced sparsity parameter is a value between 0 and 1
#' @param index a vector of integers of size p giving the group membership of 
#' each variable. By default, index=1:ncol(A) corresponds to one variable in 
#' each group 
#' @param block either 0 or 1. block==0 means that deflation is used if more 
#' than one component. A block optimisation algorithm is otherwise used 
#' that computes m components at once. By default, block=1
#' @param mu numerical vector of size m with the mu parameters (required for the block algorithms only). 
#' By default, mu_j=1/j
#' @param center a logical value indicating whether the variables should be 
#' shifted to be zero centered
#' @param scale a logical value indicating whether the variables should be 
#' scaled to have unit variance
#' @param groupsize a logical value indicating wheter the size of the groups 
#' should be taken into account. By default, groupsize=FALSE
#' @param init a matrix of size p by m to initialize the loadings matrix in 
#' the block optimisation algorithm.
#' 
#' @return \item{Z}{a p by m numerical matrix with the m group-sparse loading vectors} 
#' @return \item{Y}{a n by m numerical matrix with the m principal components}  
#' @return \item{B}{the numerical data matrix centered (if center=TRUE) and scaled 
#' (if scale=TRUE)}
#' @return \item{coef}{a numerical vector of size p with the coefficients to predict 
#' principal component scores of new observations}
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
#' center=TRUE) and scaled (if scale=TRUE) data matrix and where Z is the 
#' group-sparse loading matrix. 

#'@references 
#' M. Chavent and G. Chavent, Optimal projected variance group-sparse block PCA,
#'submitted, 2020.
#'@references
#'M. Journee, Y. Nesterov, P. Richtarik, and R. Sepulchre. Generalized power 
#'method for sparse principal component analysis. Journal of Machine Learning
#'Research, 11:517-553, 2010.
#'

#' @examples 
#' # Simulated data
#'  v1 <- c(1,1,1,1,0,0,0,0,0.9,0.9)
#'  v2 <- c(0,0,0,0,1,1,1,1,-0.3,0.3)
#'  valp <- c(200,100,50,50,6,5,4,3,2,1)
#'  A <- simuPCA(50,cbind(v1,v2),valp,seed=1)
#'  # Three group-sparse PCA algorithms
#'  index <- rep(c(1,2,3),c(4,4,2)) 
#'  Z <- groupsparsePCA(A,2,c(0.5,0.5),index,block=0)$Z #deflation
#'  Z <- groupsparsePCA(A,2,c(0.5,0.5),index)$Z # block different mu
#'  Z <- groupsparsePCA(A,2,c(0.5,0.5),index,mu=c(1,1))$Z # block same mu
#'  
#'@seealso \code{\link{sparsePCA}}, \code{\link{pev}}, \code{\link{explainedVar}}

groupsparsePCA <- function(A, m, lambda, index=1:ncol(A), block=1, mu=1/1:m,
                           groupsize=FALSE, center=TRUE, scale=TRUE,
                           init=NULL)
{
  n <- nrow(A)
  pp <- ncol(A) #number of scalar variables
  post=FALSE
  if (m > min(n-1,pp))
    stop("m must be smaller than rank(A)",call. = FALSE)
  if (length(lambda)!=m) 
    stop("lambda must be a vector of size m",call. = FALSE)
  if ((max(lambda) >1) || (max(lambda) < 0))
    stop("Values in lambda must be in [0,1]",call. = FALSE)
  if (length(index) !=ncol(A))
    stop("the length of index is not correct",call. = FALSE)
  
  norm2 <- function(x) sqrt(sum(x^2))
  polar <- function(x)
  {
    obj <- svd(x)
    obj$u %*% t(obj$v)
  }
  
  p <- length(unique(index)) #number of groups variables
  
  iter_max <- 1000   # maximum number of admissible iterations
  epsilon <- 0.0001  # accuracy of the stopping criterion
  
  
  A <- as.matrix(A)
  index <- as.factor(index)
  s <- rep(1,pp)
  
  if (center==TRUE)
  {
    A <- scale(A,center=TRUE, scale=FALSE) #center data
    g <- attr(A, "scaled:center")
  }
  if (scale==TRUE)
  {
    A <- scale(A,center=FALSE,scale=TRUE)
    s <- attr(A, "scaled:scale")*sqrt((n-1)/n)
    A <- A*sqrt(n/(n-1)) #standardized data
  }
  
  testorder <- any(unique(index)!=1:p)
  if (testorder)
  {
    A <- A[,order(index)] #column sorted from index
    rankindex <- rank(index, ties.method = "first")
    index <- sort(index)
  } 
  
  Z <- matrix(0,ncol(A),m)
  rownames(Z) <- colnames(A)  
  
  
  if ((m==1) || (m>1 && block==0)) #single-unit algorithm (deflation is used if m>1)
  {
    niter <- rep(NA, m)
    B <- A
    X <- matrix(NA,n,m)
    gamma <- rep(NA,m)
    for (comp in 1:m)            #loop on the components
    {
      ai <- lapply(split(data.frame(t(B)),index),t) # list of group variables (matrices)
      i_max <- which.max(lapply(ai,function(x){norm(x,"2")})) #norm of a matrix is the largest singular value
      gamma_max <- norm(ai[[i_max]],type="2")
      gamma[comp] <- lambda[comp]*gamma_max
      x <-svd(B)$u[,1]     #initialization point
      f <- matrix(0,iter_max,1)
      iter <- 1
      
      S <- function(v)
      {
        alpha <- sqrt(sum(v^2))
        if  (groupsize==TRUE)
          gamma[comp] <- gamma[comp]*sqrt(nrow(v))
        if (alpha <= gamma[comp])
          v[,] <- rep(0,nrow(v)) else
            v <- v*(1-gamma[comp]/alpha)
          return(v)
      }
      
      repeat
      {
        Ax <- lapply(ai,function(a){crossprod(a,x)})
        tresh <-  lapply(Ax,S)
        f[iter] <- sum(unlist(lapply(tresh,function(v){sum(v^2)})))  #cost function
        if (f[iter]==0) #the sparsity parameter is too high: all entries of the loading vector are zero
          break
        grad <-  B%*%unlist(tresh)
        x <- grad/norm2(grad)
        if (((iter>=2) && ((f[iter]-f[iter-1])/f[iter-1] < epsilon)) || (iter >= iter_max)) #stopping criterion
        {
          niter[comp] <- iter
          if (iter >= iter_max) print("Maximum number of iterations reached")
          break
        }
        iter <- iter+1
      }
      Ax <- lapply(ai,function(a){crossprod(a,x)})
      #pattern <- abs(Ax)-gamma[comp] >0
      tresh <-  lapply(Ax,S)
      z <- unlist(tresh)
      if (max(abs(z))>0)
        z <- z/norm2(z)
      #if (post==TRUE) z <- pattern_filling(A,z)
      y <- B%*%z  
      B <- B-y%*%t(z)
      Z[,comp]<-z
      X[,comp] <- y/norm2(y)
    }
  }
  if (m>1 && block==1) # block algorithm
  {
    e <- svd(A)
    d <- e$d[1:m]
    u <- e$u[,1:m]
    ai <- lapply(split(data.frame(t(A)),index),t) # list of group variables
    i_max <- which.max(lapply(ai,function(x){norm(x,"2")})) 
    gamma_max <- norm(ai[[i_max]],type="2")
    gamma <- lambda*gamma_max*d/d[1]
    
    #initialization
    if (!is.null(init)) 
    {
      x <- A %*% init %*% diag(mu ^ 2)
      x <- polar(x)
    } else
      x <- u 
    
    f <- matrix(0,iter_max,1)
    iter <- 1
    
    S <- function(v)
    {
      alpha <- apply(v,2,norm2)
      if  (groupsize==TRUE)
        gamma <- gamma * sqrt(nrow(v))
      for (comp in 1:ncol(v))
      {
        if (alpha[comp] <= gamma[comp])
          v[,comp] <- rep(0, nrow(v)) 
        else
          v[,comp] <- v[,comp]*(1-gamma[comp]/alpha[comp])
      }
      return(v)
    }
    
    repeat
    {
      Ax <- lapply(ai,function(a){crossprod(a,x)})
      tresh <-  lapply(Ax,S)
      #f[iter] <- sum(unlist(lapply(tresh,function(v){sum(v^2)})))  #cost function
      f[iter] <- sum(unlist(lapply(tresh,function(v){sum(apply(v,2,norm2)^2*mu^2)})))  #cost function
      if (f[iter]==0) #the sparsity parameter is too high: all entries of the loading vector are zero
      {
        niter <- iter
        break
      }
      for(i in seq(length(tresh))) tresh2 <- if(i == 1) tresh[[i]] else rbind(tresh2,tresh[[i]])
      grad <-  2*A%*%tresh2 %*% diag(mu^2)
      x <- polar(grad)
      if (((iter>=2) && ((f[iter]-f[iter-1])/f[iter-1] < epsilon)) || (iter >= iter_max)) #stopping criterion
      {
        niter <- iter
        if (iter >= iter_max) print("Maximum number of iterations reached")
        break
      }
      iter <- iter+1
    }
    Ax <- lapply(ai,function(a){crossprod(a,x)})
    tresh <-  lapply(Ax,S)
    for(i in seq(length(tresh))) tresh2 <- if(i == 1) tresh[[i]] else rbind(tresh2,tresh[[i]])
    z <- tresh2
    for (j in 1:m)
    {
      if (max(abs(z[,j]))>0)    
        Z[,j] <- z[,j]/norm2(z[,j])
    }
    #if (post==TRUE) Z <- pattern_filling(A,Z,mu)
    X=x
  }
  Y <- A %*% Z
  if (testorder)
  {
    Z <- Z[rankindex,] #variables back to initial order
    A <- A[,rankindex]
  }
  # coefficients for further predictions
  beta <- Z / s
  beta0 <- -t(beta) %*% matrix(g,pp,1)
  coef <- rbind(t(beta0), beta)
  colnames(coef) <- colnames(Z) <- colnames(Y) <- paste("comp",1:m, sep = "")
  res <- list(Z=Z, Y=Y, coef=coef, gamma=gamma, B=A, X=X, niter=niter)
  class(res) <- "sparsePCA"
  return(res)
}




