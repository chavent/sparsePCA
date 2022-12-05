#' @export
#' @title  Prediction of new scores
#' @description Performs scores of new observations 
#' on the principal components of sparsePCA (or groupsparsePCA). The
#' new observations must be described with the same variables than those 
#' used in sparsePCA (or groupsparsePCA).
#' @param object an object of class sparsePCA obtained with the function 
#' \code{sparsePCA} or \code{groupsparsePCA}.
#' @param newdat a numerical data matrix.
#' @param \ldots urther arguments passed to or from other methods. 
#' They are ignored in this function.
#' @return  Returns the matrix of the scores of the new observations on 
#' the principal components.
#' @seealso \code{\link{sparsePCA}}, \code{\link{groupsparsePCA}}
#' @examples 
#' #Simulated data
#' v1 <- c(1,1,1,1,0,0,0,0,0.9,0.9)
#' v2 <- c(0,0,0,0,1,1,1,1,-0.3,0.3)
#' valp <- c(200,100,50,50,6,5,4,3,2,1)
#' A <- simuPCA(50,cbind(v1,v2),valp,seed=1)
#' #Group-sparse PCA
#' index <- rep(c(1,2,3),c(4,4,2)) 
#' train <- 1:30
#' test <- 31:50
#' res <- groupsparsePCA(A[train,],2,c(0.5,0.5),index)
#' pred <- predict(res,A[test,])

predict.sparsePCA <- function(object, newdat,...)
{
  if (!inherits(object, "sparsePCA")) 
    stop("use only with \"sparsePCA\" objects")
  
  coord <- cbind(1,newdat) %*% object$coef
  colnames(coord) <- colnames(object$coef)
  rownames(coord) <- rownames(newdat)
  return(coord)			
}
