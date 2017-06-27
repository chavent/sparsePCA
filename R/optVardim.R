#' @title Explained variance of each sparse principal components
#' @description In PCA and sparse PCA, the principal components are
#' the columns of the matrix Y=BZ where B is a n times p  data matrix, and Z
#' is the p times m matrix of loadings. The matrix B can be the data is the data matrix centered and/or normalized 
#' depending on the choice of data pre-processing in the sparsePCA (or group-sparse PCA) procedure. In sparse PCA, the loadings are not necessarly orthogonal and 
#' the principal components can be correlated. In sparse PCA, the loadings are not necessarly orthogonal and 
#' the principal components can be correlated. The definition of the variance explained by each principal
#' components must then be modified. This function implements the 'optimal variance' (optVar) definition of 'explained variance'.
#' @param B a n times p  (usually centered and or scaled) data matrix. 
#' @param Z a p times m matrix of sparse loadings.
#' @return Returns the explained variance of each principal components.
#' @details The m loadings vectors in Z must be unique norm and linearly independant. 
#'  The matrix 
#' B must be the data matrix centered and/or scaled to unit variance depending on the arguments \code{center} and \code{scale} 
#' used in the function \code{sparsePCA} or \code{groupsparsePCA}. The matrix B used here must be the matrix given in output 
#' of the function used to get the matrix of sparse loadings Z (an argument of the objects of class \code{sparsePCA}).
#' @export
#' @examples 
#' # Example from Shen & Huang 2008
#'  v1 <- c(1,1,1,1,0,0,0,0,0.9,0.9)
#'  v2 <- c(0,0,0,0,1,1,1,1,-0.3,0.3)
#'  valp <- c(200,100,50,50,6,5,4,3,2,1)
#'  A <- simuPCA(50,cbind(v1,v2),valp,seed=1)
#'  Z <- sparsePCA(A,2,c(0.5,0.5))$Z #deflation algo
#'  B <- sparsePCA(A,2,c(0.5,0.5))$B
#'  optVardim(B,Z)
#'@references 
#'\itemize{
#'\item M. Chavent and G. Chavent, Group-sparse block PCA and explained variance, arXiv:1705.00461
#'}
#'@seealso \code{\link{sparsePCA}}, \code{\link{groupsparsePCA}},\code{\link{pev}}
#'
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


