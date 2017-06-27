#' @title proportion of correct or incorrect zero
#' @description Returns the proportion of zero correctly or incorrectly found in each colum of the loadings matrix Z.
#' @param Ztrue the true sparse loadings matrix.
#' @param Z the obtained sparse loading matrix.
#' @param method either "correct" or "incorrect".
#' @details The proportion of zero correctly found is the true positive rate (tpr) i.e. the proportion of true zero
#' put to zero in Z. The proportion of zero incorrectly found is the false positive rate (fpr) i.e. the proportion
#' of zero in Z which are not true zero.
#' @examples 
#' # Example from Shen & Huang 2008
#'  v1 <- c(1,1,1,1,0,0,0,0,0.9,0.9)
#'  v2 <- c(0,0,0,0,1,1,1,1,-0.3,0.3)
#'  Ztrue <- cbind(v1,v2)
#'  valp <- c(200,100,50,50,6,5,4,3,2,1)
#'  n <- 50
#'  A <- simuPCA(n,Ztrue,valp,seed=1)
#'  Z <- sparsePCA(A,2,c(0.5,0.5),block=1)$Z
#'  truth(Ztrue,Z) #tpr
#'  truth(Ztrue,Z,method="incorrect") #fpr
#' @export
truth <- function(Ztrue,Z,method="correct")
{
  if (!(method %in% c("correct","incorrect")))
    stop("the names of the method is not correct", call. = FALSE)
  TPR <- function(y_pred, y) {
    VP = sum((y_pred == TRUE) & (y == TRUE))
    return(VP/sum(y == TRUE))
  }
  
  TNR <- function(y_pred, y) {
    VN = sum((y_pred == FALSE) & (y == FALSE))
    return(VN/sum(y == FALSE))
  }
  P1 <- !Ztrue #true positions of 0
  P2 <- !Z #predicted positions of 0 
  pourc <- rep(NA,ncol(Ztrue))
  for (j in 1:ncol(P1))
  {
    if (method=="correct")
      pourc[j] <- TPR(P2[,j],P1[,j])
    if (method=="incorrect")
      pourc[j] <- 1-TNR(P2[,j],P1[,j])
  }
  return(pourc)
}


#' @title RV coefficient
#' @description The RV coefficient is a measures the proximity between two matrices (subspaces spanned by the columns ofthe matrices). This measure
#' is normalized and takes its values between 0 and 1.
#' @param X a matrix of dimension n times p
#' @param Y a matrix of dimension n times q
#' @export
RVcoef <- function(X,Y)
{
  sum((t(X)%*%Y)^2)/(sqrt(sum((t(X)%*%X)^2))*sqrt(sum((t(Y)%*%Y)^2)))
}

#' @title  Volume coefficient
#' @description The volume of the parallepipede constructed on the columns of the matrix Y=BZ. This volume measures the orthogonality of the columns of Y.
#' It is normalized and takes its values between 0 and 1.
#' @param B a n times p (centered or standardized) data matrix.
#' @param Z a p times m matrix of sparse loadings.
#' @export
Vol <- function(B,Z)
{
  norm2 <- function(x) sqrt(sum(x^2))
  if (sum(abs(Z))>0)
  {
    sel <- apply(abs(Z),2,sum) > 0
    Y <- B%*%Z[,sel,drop=FALSE]
    nor <- apply(Y,2,norm2)
    Y <- Y[,order(nor,decreasing=TRUE)]
    val <- abs(prod(diag(qr.R(qr(Y)))))/prod(nor)
  } else val <- NA
  return(val)
    
}

#Degree of sparsity
degsp <- function(Z)
{
  res <- apply(Z==0,2,sum)
  return(res)
}

#Frobenius norm
frobenius <- function(B)
{  
  obj <- svd(B)
  sum(obj$d^2)
}


#polar matrix U
polar <- function(x)
{
  obj <- svd(x)
  obj$u%*%t(obj$v)
}