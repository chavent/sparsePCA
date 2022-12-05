#' @title Variance explained by non orthogonal principal components
#' @export
#' 
#' @description  The loadings in sparse and group-sparse PCA
#' are not necessarly orthogonal and new definitions of explained variance must
#' be used. This function implements six different definitions of explained 
#' variance:
#' the subspace variance (subspVar), the adjusted variance (adjVar), 
#' the optimal projected variance (optVar)
#' the polar variance (polVar), the 
#' QR normalized variance (QRnormVar) and the UP normalized variance (UPnormVar).
#' 
#' @param B a numerical data matrix  (usually centered and/or scaled, see details).
#' @param Z a numerical loading matrix.
#' @param method the name of the explained variance: "subspVar", "optVar", 
#' "polVar", "adjVar","QRnormVar","UPnormVar".
#' 
#' @return Returns the variance explained by the components.
#' 
#' @details 
#' The principal components are defined by
#' Y=BZ where B is the centered and/or scaled data matrix and Z is the sparse 
#' loading matrix. The argument \code{B} must then be consistent with the 
#' pre-processing step in \code{\link{sparsePCA}} and \code{\link{groupsparsePCA}}.
#' The loading vectors in \code{Z} must be unique norm and linearly independant. 
#' 
#' The "adjVar" method is first introduced in Zou et al. (2006), the "subspVar" 
#' method is introduced in Shen & Huang (2008) and the "optVar", "QRnormVar" 
#' and "UPnormVar" methods are introduced in Chavent & Chavent (2020).
#' The six variance definitions are compared in Chavent & Chavent (2020)
#' 

#'  
#'@references M. Chavent and G. Chavent, Optimal projected variance group-sparse
#' block PCA, submitted, 2020.
#' @references  H. Zou, T. Hastie, R. Tibshirani. Sparse principal component 
#' analysis. Journal of computational and graphical statistics, 15(2):265-286, 2006.
#' @references H. Shen and J.Z. Huang. Sparse principal component analysis via 
#' regularized low rank matrix approximation. Journal of multivariate analysis, 
#' 99(6):1015-1034, 2008.
#' 
#'@seealso \code{\link{pev}}, \code{\link{sparsePCA}}, \code{\link{groupsparsePCA}} 

#' @examples 
#' # Example from Shen & Huang 2008
#'  v1 <- c(1,1,1,1,0,0,0,0,0.9,0.9)
#'  v2 <- c(0,0,0,0,1,1,1,1,-0.3,0.3)
#'  valp <- c(200,100,50,50,6,5,4,3,2,1)
#'  A <- simuPCA(50,cbind(v1,v2),valp,seed=1)
#'  Z <- sparsePCA(A,2,c(0.5,0.5))$Z #deflation
#'  B <- sparsePCA(A,2,c(0.5,0.5))$B
#'  namevar <- c("subspVar", "optVar", "polVar", "adjVar","QRnormVar","UPnormVar")
#'  vectvar <- rep(NA,length(namevar))
#'  names(vectvar) <- namevar
#'  for (j in 1:length(namevar)) vectvar[j] <- explainedVar(B,Z,method=namevar[j])
#'  
explainedVar <- function(B, Z, method="optVar")
{
  if (!(method %in% c("subspVar", "optVar", "polVar", "adjVar", "QRnormVar", "UPnormVar")))
    stop("the names of the method is not correct", call. = FALSE)
  
  if (nrow(Z) != ncol(B))
    stop("the dimensions of B and/or Z are not correct", call. = FALSE)
  
  if (!(is.matrix(B) && is.matrix(Z)))
    stop("B and Z must be objects of class matrix", call. = FALSE)
  
  Tr <- function(x) sum(diag(x))
  norm2 <- function(x) sqrt(sum(x^2))
  polar <- function(x)
  {
    obj <- svd(x)
    obj$u %*% t(obj$v)
  }
  
  if (sum(abs(Z))>0) #at least one non zero loading
  {
    sel <- apply(abs(Z), 2, sum) > 0
    Y <- B%*%Z[, sel, drop=FALSE]
    
    if (method=="subspVar")
    {
      S <- chol2inv(chol(crossprod(Z[,sel]))) # (t(Z)%*%Z)^{-1}
      val <- Tr(crossprod(Y)%*%S) #Trace(t(Y)%*%Y%*%S)
    }
    
    if (method=="adjVar")
    {
      Y <- Y[, order(apply(Y,2,norm2), decreasing=TRUE)]
      R <- qr.R(qr(Y)) # QR decomposition of Y
      val = sum(diag(R)^2) #Trace(diag(R)^2)
    }
    
    if (method=="polVar")
    {
      if (sum(sel !=0)==1) val <- t(Y)%*%Y
      else {
        #C <- t(Y) %*% Y
        #e <- eigen(C,symmetric=TRUE)
        #V <- e$vectors
        #srC <- V %*% diag(sqrt(e$values)) %*% t(V) # (t(Y) %*% Y)^{1/2}
        X <- polar(Y)
        M <- t(X) %*% Y
        val <- sum(diag(M)^2)
      }
    }
    
    iter <- NULL
    if (method=="optVar")
    {
      if (sum(sel !=0)==1) val <- t(Y)%*%Y
      else {
        X0 <- polar(Y)
        iter <- 0
        repeat
        {
          X <- polar(Y %*% diag(diag(t(X0)%*%Y)))
          iter <- iter+1
          if (isTRUE(all.equal(X0 , X)))
            break else X0=X
        }
        M <- t(X) %*% Y
        val <- sum(diag(M)^2)
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


