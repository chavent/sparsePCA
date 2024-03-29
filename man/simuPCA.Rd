% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simuPCA.R
\name{simuPCA}
\alias{simuPCA}
\title{Simulation from a PCA model}
\usage{
simuPCA(n = 30, Ztrue, valp, seed = FALSE)
}
\arguments{
\item{n}{the number of observations}

\item{Ztrue}{a numerical matrix of size p by m with the underlying loading 
vectors (m first eigenvectors of the underlying covariance
matrix)}

\item{valp}{a numerical vector of size p with the eigenvalues of 
the underlying covariance matrix}

\item{seed}{a seed value to replicate the simulation. By default, seed=FALSE}
}
\value{
Returns a numerical data matrix of size n by p
}
\description{
This function  simulates data from a PCA model. 
The data matrix is drawn from a centered multivariate normal 
distribution where the covariance matrix is constructed 
knowing its p eigenvalues and its m < p first eigenvectors 
(see Chavent & Chavent 2020).
}
\examples{
# Example from Journee & al. 2010
 v1 <- c(rep(1/sqrt(10),10),rep(0,490))
 v2 <- c(rep(0,10), rep(1/sqrt(10),10),rep(0,480))
 valp <- c(400,300,rep(1,498))
 n <- 100
 A <- simuPCA(n,cbind(v1,v2),valp,seed=1)
 svd(A)$d[1:3]^2 #eigenvalues of the empirical covariance matrix.
 
# Example from Shen & Huang 2008
 v1 <- c(1,1,1,1,0,0,0,0,0.9,0.9)
 v2 <- c(0,0,0,0,1,1,1,1,-0.3,0.3)
 valp <- c(200,100,50,50,6,5,4,3,2,1)
 n <- 100
 A <- simuPCA(n,cbind(v1,v2),valp,seed=1)
 svd(A)$d^2 #eigenvalues of the empirical covariance matrix (times n).
 
# Example from Chavent & Chavent 2020
 valp <- c(c(200,180,150,130),rep(1,16))
 data("Ztrue")
 n <- 100
 A <-simuPCA(n,Ztrue,valp,seed=1)
 svd(A)$d^2 #eigenvalues of the empirical covariance matrix.
 
}
\references{
M. Chavent and G. Chavent, Optimal projected variance group-sparse block PCA, 
submitted, 2020.

M. Journee, Y. Nesterov, P. Richtarik, and R. Sepulchre. Generalized power 
method for sparse principal component analysis. Journal of Machine Learning
Research, 11:517-553, 2010.

H. Shen and J.Z. Huang. Sparse principal component analysis via 
regularized low rank matrix approximation. Journal of multivariate analysis, 
99(6):1015-1034, 2008.
}
\seealso{
\code{\link{sparsePCA}}, \code{\link{groupsparsePCA}}
}
