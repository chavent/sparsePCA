% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/explainedVar.R
\name{explainedVar}
\alias{explainedVar}
\title{Total explained variance of sparse principal components}
\usage{
explainedVar(B, Z, method = "optVar")
}
\arguments{
\item{B}{a n times p  (usually centered and or scaled) data matrix.}

\item{Z}{a p times m matrix of sparse loadings.}

\item{method}{the name of the explained variance: "subspVar","adjVar","optVar","QRnormVar","UPnormVar".}
}
\value{
Returns the explained variance of the principal components.
}
\description{
In PCA and sparse PCA, the principal components are
the columns of the matrix Y=BZ where B is a n times p  data matrix, and Z
is the p times m matrix of loadings. The matrix B can be the data is the data matrix centered and/or normalized 
depending on the choice of data pre-processing in the sparsePCA (or group-sparse PCA) procedure. In sparse PCA, the loadings are not necessarly orthogonal and 
the principal components can be correlated. The definition of the variance explained by the principal
components must then be modified. This function implements five different definitions of 'explained variance':
the subspace variance (subspVar), the adjusted variance (adjVar), the optimal variance (optVar), 
the QR normalized variance (QRnormVar) and the UP normalized variance (UPnormVar).
}
\details{
The m loadings vectors in Z must be unique norm and linearly independant.  The matrix 
B must be the data matrix centered and/or scaled to unit variance depending on the arguments \code{center} and \code{scale} 
used in the function \code{sparsePCA} or \code{groupsparsePCA}. The matrix B used here must be the matrix given in output 
of the function used to get the matrix of sparse loadings Z (an argument of the objects of class \code{sparsePCA}).

If method="subspVar", the explained variance is the variance of the projection of B 
on the subspace spanned by the m loading vectors. If has been introduced in Shen & Huand (2008).

If method="adjVar", the adjusted explained variance is that introduced by Zou et al. (2006).
}
\examples{
# Example from Shen & Huang 2008
 v1 <- c(1,1,1,1,0,0,0,0,0.9,0.9)
 v2 <- c(0,0,0,0,1,1,1,1,-0.3,0.3)
 valp <- c(200,100,50,50,6,5,4,3,2,1)
 A <- simuPCA(50,cbind(v1,v2),valp,seed=1)
 Z <- sparsePCA(A,2,c(0.5,0.5))$Z #deflation algo
 B <- sparsePCA(A,2,c(0.5,0.5))$B
 namevar <- c("subspVar","adjVar","optVar","QRnormVar","UPnormVar")
 vectvar <- rep(NA,length(namevar))
 names(vectvar) <- namevar
 for (j in 1:length(namevar)) vectvar[j] <- explainedVar(B,Z,method=namevar[j])
}
\references{
\itemize{
\item M. Chavent and G. Chavent, Group-sparse block PCA and explained variance, arXiv:1705.00461
}
}
\seealso{
\code{\link{sparsePCA}}, \code{\link{groupsparsePCA}},\code{\link{pev}}
}

