#' @title Total explained variance of sparse principal components
#' @description In PCA and sparse PCA, the principal components are
#' the columns of the matrix Y=BZ where B is a n times p centered (or standardized) data matrix, and Z
#' is the p times m matrix of loadings. In sparse PCA, the loadings are not necessarly orthogonal and 
#' the principal components can be correlated. The definition of the variance explained by the principal
#' components must then be modified. This function implements five different definitions of 'explained variance':
#' the subspace variance (subspVar), the adjusted variance (adjVar), the optimal variance (optVar), 
#' the QR normalized variance (QRnormVar) and the UP normalized variance (UPnormVar).
#' @param B a n times p  (centered or standardized) data matrix.
#' @param Z a p times m matrix of sparse loadings.
#' @param method the name of the explained variance: "subspVar","adjVar","optVar","QRnormVar","UPnormVar".
#' @return Returns the explained variance of the principal components.
#' @details The m loadings vectors in Z must be unique norm and linearly independant. 
#' 
#' If method="subspVar", the explained variance is the variance of the projection of B 
#' on the subspace spanned by the m loading vectors. If has been introduced in Shen & Huand (2008).
#' 
#' If method="adjVar", the adjusted explained variance is that introduced by Zou et al. (2006).
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
#'  namevar <- c("subspVar","adjVar","optVar","QRnormVar","UPnormVar")
#'  vectvar <- rep(NA,length(namevar))
#'  names(vectvar) <- namevar
#'  for (j in 1:length(namevar)) vectvar[j] <- explainedVar(B,Z,method=namevar[j])

explainedVar <- function(B,Z,method="optVar")
{
  if (!(method %in% c("subspVar","adjVar","optVar","QRnormVar","UPnormVar")))
    stop("the names of the method is not correct", call. = FALSE)
  if (nrow(Z) != ncol(B))
    stop("the dimension of B and/or Z are not correct",call. = FALSE)
  if (!(is.matrix(B) && is.matrix(Z)))
    stop("B and Z must be objects of class matrix",call. = FALSE)
  Tr <- function(x) sum(diag(x))
  norm2 <- function(x) sqrt(sum(x^2))
  
  if (sum(abs(Z))>0) #at least one non zero loading
  {
    sel <- apply(abs(Z),2,sum) > 0
    Y <- B%*%Z[,sel,drop=FALSE]
    if (method=="subspVar")
    {
      S <- chol2inv(chol(crossprod(Z[,sel]))) # (t(Z)%*%Z)^{-1}
      val <- Tr(crossprod(Y)%*%S) #Trace(t(Y)%*%Y%*%S)
    }
    if (method=="adjVar")
    {
      Y <- Y[,order(apply(Y,2,norm2),decreasing=TRUE)]
      R <- qr.R(qr(Y)) # QR decomposition of Y
      val=sum(diag(R)^2) #Trace(R^2)
    }
    if (method=="optVar")
    {
      if (sum(sel !=0)==1) val <- t(Y)%*%Y
      else {
        C <- t(Y) %*% Y
        e <- eigen(C,symmetric=TRUE)
        V <- e$vectors
        srC <- V %*% diag(sqrt(e$values)) %*% t(V) # (t(Y) %*% Y)^{1/2}
        val <- sum(diag(srC)^2)
      }
    }
    if (method=="QRnormVar")
    {
      Y <- Y[,order(apply(Y,2,norm2),decreasing=TRUE)]
      R <- qr.R(qr(Y)) # QR decomposition of Y
      TT <- Z[,sel,drop=FALSE] %*% solve(R) # TT=ZR^{-1}
      val <- sum(1/diag(t(TT)%*%TT)) 
    }
    if (method=="UPnormVar")
    {
      C <- t(Y) %*% Y
      e <- eigen(C,symmetric=TRUE)
      V <- e$vectors
      if (length(e$values)==1)
        TT <- Z[,sel,drop=FALSE] /sqrt(e$values)
      else {
        invsrC <- V %*% diag(1/sqrt(e$values)) %*% t(V) # (t(Y) %*% Y)^{-1/2}
        TT <- Z[,sel,drop=FALSE] %*% invsrC 
      }
      val <- sum(1/diag(t(TT)%*%TT)) 
    }
  } else val <- 0
  return(val)
}





