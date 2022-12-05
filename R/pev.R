#' @title Proportion of variance explained by non orthogonal principal components
#' @export
#' 
#' @description This function calculates the proportion of variance explained by
#' not necessarily orthogonal principal components, using the optimal projected
#' variance (optVar). 
#' 
#' @param x a object of class \code{sparsePCA}.
#' @return The proportion of variance explained by each principal components. 
#' 
#' @details 
#' The object \code{x} is the output of \code{\link{sparsePCA}} or 
#' \code{\link{groupsparsePCA}}. The  variance of the not necessarily 
#' orthogonal principal components in \code{x}
#' is calculated with the optimal projected variance definition (see Chavent and Chavent, 2020)
#' and is divided by the total variance of the data to get a proportion of
#' explained variance (pev).
#' 
#'@references M. Chavent and G. Chavent, Optimal projected variance group-sparse
#' block PCA, submitted, 2020.
#' 
#' @seealso   \code{\link{optVardim}}, \code{\link{explainedVar}}, 
#' \code{\link{sparsePCA}}, \code{\link{groupsparsePCA}}, 

#' @examples 
#' data(protein)
#' x <- sparsePCA(protein, 2, c(0.5,0.5))
#' pev(x)
#' x <- groupsparsePCA(protein, 2, index=1:ncol(protein), c(0.5,0.5))
#'pev(x)


pev <- function(x)
  {
  if (!inherits(x, "sparsePCA")) 
      	stop("use only with \"sparsePCA\" objects")
  
  return(optVardim(x$B,x$Z)/sum(svd(x$B)$d^2))
}