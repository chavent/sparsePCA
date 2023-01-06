#' @export
#' @name sparsePCAmix
#' @title Sparse principal component analysis of mixed data 
#' @description Performs sparse principal component analysis of  a set of 
#' individuals (observations) described by a mixture of qualitative and 
#' quantitative variables. sparsePCAmix includes ordinary sparse principal component 
#' analysis (PCA) and sparse multiple correspondence analysis (MCA) as special cases.
#' @param X.quanti a numeric matrix of data.
#' @param X.quali a categorical matrix of data.
#' @param m number of sparse components.
#' @param lambda a vector of dimension m with reduced sparsity parameters (in relative value with respect to the theoretical upper bound). 
#' Each reduced sparsity parameter is a value between 0 and 1.
#' @param block either 0 or 1. block==0 means that deflation is used if more than one component need to be computed.
#' A block algorithm is otherwise used, that computes m components at once. By default, block=1.
#' @param mu vector of dimension m with the mu parameters (required for the block algorithms only). 
#' By default, mu_j=1/j
#' @param groupsize a logical value indicating wheter the size of the groups should be taken into account.
#' @param rename.level boolean, if TRUE all the levels of the qualitative 
#' variables are renamed as follows: "variable_name=level_name". 
#' This prevents to have identical names of the levels.
#' @return \item{V}{the p times m matrix that contains the m sparse loading vectors} 
#' @return \item{scores}{the n times m matrix that contains the m principal components}  
#' @return \item{pev}{the proportion of variance (calculated with 'optVar') explained by the components.}
#' @return \item{varsel}{a list with the name of the variables selected in each dimension.}
#' @return \item{degsp}{the degree of sparsity of each component (number of selected variables).}
#' @return \item{Z}{the pre-processed data matrix}
#' @return \item{w}{the vector of the weights of the columns of Z.}
#' @details The pre-processed data matrix Z is X.quanti standardized (centered and 
#' reduced by standard deviations) concatenated with the indicator matrix of X.quali 
#' centered. The principal components (the scores) are given by the matrix F=ZMV where M=diag(w) 
#' is the diagonal metric of the weights of the columns of Z.  
#' In sparse PCA, the loadings are not necessarly orthogonal and 
#' the principal components can be correlated. The definition of the variance explained 
#' by each principal components must then be modified. Here the pev (proportion of explained
#' variance) of the PCs is calculated with the 'optimal variance' (optVar) definition 
#' of 'explained variance'.
#'@references 
#'\itemize{
#'\item M. Chavent and G. Chavent, Group-sparse block PCA and explained variance, arXiv:1705.00461
#'\item M. Chavent, V. Kuentz-Simonet, A. Labenne, J. Saracco. Multivariate Analysis of Mixed Data: The R Package PCAmixdata. Electronic Journal of Applied Statistical Analysis. ⟨hal-01662595⟩
#'}

sparsePCAmix<- function (X.quanti = NULL, X.quali = NULL, m = 2, lambda,
                         block = 1, mu=1/1:m, groupsize=FALSE, rename.level=FALSE)
{
  cl <- match.call()
  rec <- PCAmixdata::recod(X.quanti, X.quali, rename.level)
  n <- rec$n
  p <- rec$p
  p1 <- rec$p1
  p2 <- p - p1

  #Build the metrics 
  N <- rep(1/n, n)
  M1 <- M2 <- NULL 
  if (p1!=0)  M1 <- rep(1,p1) 
  if (p2!=0){
    ns <- apply(rec$G, 2, sum)
    M2 <- n/ns
   }
  M <- c(M1,M2) # weights of the columns
  names(M) <- colnames(rec$W)
  
  #Recoding step
  Z <- rec$W # first recoding step of PCAmix
  Ztilde <- sweep(Z, 2, sqrt(M), FUN = "*") 
  Ztilde <- sweep(Ztilde, 1, sqrt(N), FUN = "*")

  #Group-sparse PCA
  index <- rec$indexj
  res <- sparsePCA::groupsparsePCA(Ztilde, index, lambda=lambda , 
                                   m=m,block=block,
                                   groupsize=groupsize, center=TRUE, 
                                   scale=FALSE)
  Vtilde <- res$Z
  V <- sweep(Vtilde, 1, sqrt(M), FUN = "/")
  
  #Principal components
  FF <- Z%*%diag(M)%*%V
  
  colnames(FF) <- colnames(V) <-paste("dim", 1:m, sep = " ")
  
  #Explained variance of the PCs
  pev <- sparsePCA::pev(res)
  degsp <- apply(V!=0, 2,function(x) (length(unique(rec$indexj[x]))))
  varsel <- lapply(data.frame(V!=0), function(x) (colnames(rec$X)[unique(rec$indexj[x])]))
  
  res <- list(call = cl, 
              V = V, 
              scores = FF,
              pev=pev,
              varsel=varsel,
              degsp=degsp,
              w=M,
              Z = Z,
              Vtilde=Vtilde
  )
  return(res)
}

optVardim_metrics <- function(B,V,N,M)
{
  val <- rep(0,ncol(V))
  if (sum(abs(V))>0)
  {
    sel <- apply(abs(V),2,sum) > 0
    Y <- B%*%diag(M)%*%V[,sel,drop=FALSE]
    #Y <- B%*%Z[,,drop=FALSE]
    if (sum(sel !=0)==1) val[sel] <- t(Y)%*%Y
    else {
    C <- t(Y) %*% diag(N)%*%Y
    e <- eigen(C,symmetric=TRUE)
    v <- e$vectors
    srC <- v %*% diag(sqrt(e$values)) %*% t(v) # (t(Y) %*% Y)^{1/2}
    val[sel] <- diag(srC)^2
    }
  } else val <- rep(0,ncol(V))
  return(val)
}
