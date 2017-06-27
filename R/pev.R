#' @title Proportion of explained variance of sparse principal components
#' @description The proportion of explained variance of the m principal components i.e. the optimal variance (optVar) divided by
#' the total variance of the data (the Frobenius norm of A).
#' @param x a object of class \code{sparsePCA}.
#' @return \item{pev}{the proportion of explained variance of the m principal components i.e. the optimal variance (optVar) divided 
#' the total variance of the data (the Frobenius norm of A).} 
#' @references 
#'\itemize{
#'\item M. Chavent & G. Chavent, Group-sparse block PCA and explained variance, arXiv:1705.00461
#'}
#' @export
#' @examples 
#' data(protein)
#' x <- sparsePCA(protein,2,c(0.5,0.5))
#' pev(x)
#' x <- groupsparsePCA(protein,2,index=1:ncol(protein),c(0.5,0.5))
#'pev(x)
#'@seealso \code{\link{sparsePCA}}, \code{\link{groupsparsePCA}}, \code{\link{explainedVar}}, \code{\link{optVardim}}, 
pev <- function(x)
  {
  if (!inherits(x, "sparsePCA")) 
      	stop("use only with \"sparsePCA\" objects")
  return(optVardim(x$B,x$Z)/sum(svd(x$B)$d^2))
}