#' @title Explained variance of sparse principal components
#' @description In PCA and sparse PCA, the principal components are
#' the columns of the matrix Y=BZ where B is a n times p centered (or standardized) data matrix, and Z
#' is the p times m matrix of loadings. In sparse PCA, the loadings are not necessarly orthogonal and 
#' the principal components can be correlated. The definition of the variance explained by each principal
#' components must then be modified. This function implements the 'optimal variance' (optVar) definition of 'explained variance'.
#' @param B a n times p  (centered or standardized) data matrix.
#' @param Z a p times m matrix of sparse loadings.
#' @return Returns the explained variance of each principal components.
#' @details The m loadings vectors in Z must be unique norm and linearly independant. 
#' 
#' @export
#' @examples 
#' # Example from Shen & Huang 2008
#'  v1 <- c(1,1,1,1,0,0,0,0,0.9,0.9)
#'  v2 <- c(0,0,0,0,1,1,1,1,-0.3,0.3)
#'  valp <- c(200,100,50,50,6,5,4,3,2,1)
#'  A <- simusparsePCA(50,cbind(v1,v2),valp,seed=1)
#'  Z <- sparsePCA(A,2,c(0.5,0.5))$Z #deflation algo
#'  B <- sparsePCA(A,2,c(0.5,0.5))$B
#'  optVardim(B,Z)

optVardim <- function(B,Z)
{
  val <- rep(0,ncol(Z))
  if (sum(abs(Z))>0)
  {
    sel <- apply(abs(Z),2,sum) > 0
    Y <- B%*%Z[,sel,drop=FALSE]
    #Y <- B%*%Z[,,drop=FALSE]
    if (sum(sel !=0)==1) val[sel] <- t(Y)%*%Y
    else {
    C <- t(Y) %*% Y
    e <- eigen(C,symmetric=TRUE)
    V <- e$vectors
    srC <- V %*% diag(sqrt(e$values)) %*% t(V) # (t(Y) %*% Y)^{1/2}
    val[sel] <- diag(srC)^2
    }
  } else val <- rep(0,ncol(Z))
  return(val)
}


