#' @title Optimal projected variance of non orthogonal principal components
#' @export
#' 
#' @description This function calculates the variance explained by
#' not necessarily orthogonal principal components, using the optimal projected 
#' variance (optVar). 
#' 
#' @param B a numerical data matrix  (usually centered and/or scaled, see details).
#' @param Z a numerical loading matrix.
#' 
#' @return Returns the optimal projected variance of each principal components.
#' @details 
#' The principal components are defined by
#' Y=BZ where B is the centered and/or scaled data matrix and Z is the sparse 
#' loading matrix. The argument \code{B} must then be consistent with the 
#' pre-processing step in \code{\link{sparsePCA}} and \code{\link{groupsparsePCA}}.
#' The loading vectors in \code{Z} must be unique norm and linearly independant.
#' 

#' @examples 
#' data(protein)
#' B <- sparsePCA(protein, 2, c(0.5,0.5))$B
#' Z <- sparsePCA(protein, 2, c(0.5,0.5))$Z
#' optVardim(B,Z)
#'  
#'@references M. Chavent and G. Chavent, Optimal projected variance group-sparse
#' block PCA, submitted, 2020.
#' 
#'@seealso \code{\link{pev}}, \code{\link{explainedVar}}, 
#'\code{\link{sparsePCA}}, \code{\link{groupsparsePCA}}
#'
optVardim <- function(B,Z)
{
  val <- rep(0,ncol(Z))
  if (sum(abs(Z))>0)
  {
    sel <- apply(abs(Z),2,sum) > 0
    Y <- B%*%Z[,sel,drop=FALSE]
    if (sum(sel !=0)==1) val[sel] <- t(Y)%*%Y
    else {
      X0 <- polar(Y)
      repeat
      {
        X <- polar(Y %*% diag(diag(t(X0)%*%Y)))
        if (isTRUE(all.equal(X0 , X)))
          break else X0=X
      }
      M <- t(X) %*% Y
      val[sel] <- diag(M)^2
    }
  } else val <- rep(0,ncol(Z))
  return(val)
}

