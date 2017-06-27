#' @title Group-sparse PCA 
#' @description Group-sparse principal component analysis (PCA). Variables are arganized in groups and loadings of variables within the same
#' group are set to zero simultaneously. Three algorithms are proposed for the determination of the group-sparse loadings : on based on a deflation procedure and two based on block optimisation.
#' @param A a n times p  numerical data matrix.
#' @param m number of components.
#' @param lambda a vector of dimension m with reduced sparsity parameters (in relative value with respect to the theoretical upper bound). 
#' Each reduced sparsity parameter is a value between 0 and 1.
#' @param index vector  which  defines  the  grouping  of  the  variables. Components sharing 
#' the same number build a group. By default, index=1:ncol(A) corresponds to one variable in each group i.e. (non group) sparse PCA. 
#' @param block either 0 or 1. block==0 means that deflation is used if more than one component need to be computed.
#' A block algorithm is otherwise used, that computes m components at once. By default, block=1.
#' @param mu vector of dimension m with the mu parameters (required for the block algorithms only). 
#' By default, mu_j=1/j
#' @param post If TRUE the loadings are obtained with post-processing.
#' @param center a logical value indicating whether the variables should be shifted to be zero centered.
#' @param scale a logical value indicating whether the variables should be scaled to have unit variance.
#' @param init a matrix of dimension p times m to initialize the loadings matrix Z in the algorithm.
#' @return \item{Z}{the p times m matrix that contains the m sparse loading vectors} 
#' @return \item{Y}{the n times m matrix that contains the m principal components}  
#' @return \item{B}{the data matrix centered (if center=TRUE) and scaled (if scale=TRUE)}
#' @details The principal components are given by the matrix Y=BZ  where B is the matrix A which has been centered (if scale=FALSE) or standardized (if scale=TRUE) and Z is the matrix of the sparse loading vectors. 
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
#'  # 10 scalar variables and 3 group variables of size 4,4,2
#'  index <- rep(c(1,2,3),c(4,4,2)) 
#'  Z <- groupsparsePCA(A,2,c(0.5,0.5),index,block=0)$Z #deflation
#'  Z <- groupsparsePCA(A,2,c(0.5,0.5),index)$Z # block different mu
#'  Z <- groupsparsePCA(A,2,c(0.5,0.5),index,mu=c(1,1))$Z # block same mu
#'@seealso \code{\link{sparsePCA}}, \code{\link{pev}}

groupsparsePCA <- function(A, m, lambda, index=1:ncol(A), block=1, mu=1/1:m,post=FALSE,center=TRUE,scale=TRUE,init=NULL)
{
  n <- nrow(A)
  pp <- ncol(A) #number of scalar variables
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
  
  Z <- matrix(0,ncol(A),m)
  rownames(Z) <- colnames(A)
  A <- as.matrix(A)
  index <- as.factor(index)

  if (center==TRUE)
    A <- scale(A,center=TRUE, scale=FALSE) #center data
  if (scale==TRUE)
    A <- scale(A,center=FALSE,scale=TRUE)*sqrt(n/(n-1)) #scale data
 
  if ((m==1) || (m>1 && block==0)) #single-unit algorithm (deflation is used if m>1)
  {
    B <- A
    X <- matrix(NA,n,m)
    for (comp in 1:m)            #loop on the components
    {
      ai <- lapply(split(data.frame(t(B)),index),t) # list of group variables (matrices)
      i_max <- which.max(lapply(ai,function(x){norm(x,"2")})) #norm of a matrix is the largest singular value
      gamma_max <- norm(ai[[i_max]],type="2")
      gamma <- lambda[comp]*gamma_max
      x <-svd(B)$u[,1]     #initialization point
      f <- matrix(0,iter_max,1)
      iter <- 1
      S <- function(v)
      {
        alpha=sqrt(sum(v^2))
        if (alpha > 0) v/alpha*((alpha>=gamma)*(alpha-gamma))
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
          if (iter >= iter_max) print("Maximum number of iterations reached")
          break
        }
        iter <- iter+1
      }
      Ax <- lapply(ai,function(a){crossprod(a,x)})
      #pattern <- abs(Ax)-gamma >0
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
      sweep(v,2,alpha,"/")%*%diag((alpha >= gamma)*(alpha-gamma))
    }
    repeat
    {
      Ax <- lapply(ai,function(a){crossprod(a,x)})
      tresh <-  lapply(Ax,S)
      #f[iter] <- sum(unlist(lapply(tresh,function(v){sum(v^2)})))  #cost function
      f[iter] <- sum(unlist(lapply(tresh,function(v){sum(apply(v,2,norm2)^2*mu^2)})))  #cost function
      if (f[iter]==0) #the sparsity parameter is too high: all entries of the loading vector are zero
        break
      for(i in seq(length(tresh))) tresh2 <- if(i == 1) tresh[[i]] else rbind(tresh2,tresh[[i]])
      grad <-  2*A%*%tresh2 %*% diag(mu^2)
      x <- polar(grad)
      if (((iter>=2) && ((f[iter]-f[iter-1])/f[iter-1] < epsilon)) || (iter >= iter_max)) #stopping criterion
      {
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
  res <- list(Z=Z,Y=Y,gamma=gamma,B=A,X=X)
  class(res) <- "sparsePCA"
  return(res)
}




