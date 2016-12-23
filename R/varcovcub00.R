#' @title Variance-covariance matrix of a CUB model without covariates
#' @description Compute the variance-covariance matrix of parameter estimates of a CUB model without covariates.
#' @aliases varcovcub00
#' @usage varcovcub00(m, ordinal, pai, csi)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param pai Uncertainty parameter
#' @param csi Feeling parameter
#' @export varcovcub00
#' @details The function checks if the variance-covariance matrix is positive-definite: if not, 
#' it returns a warning message and produces a matrix with NA entries.
#' @seealso \code{\link{probcub00}}
#' @keywords internal
#' @references
#' Piccolo D. (2006), Observed Information Matrix for MUB Models. \emph{Quaderni di Statistica},
#'  \bold{8}, 33--78,
#' @examples
#' data(univer)
#' m<-7
#' ordinal<-univer[,12]
#' pai<-0.87
#' csi<-0.17
#' varmat<-varcovcub00(m, ordinal, pai, csi)


varcovcub00 <-
function(m,ordinal,pai,csi){
  vvi<-(m-ordinal)/csi-(ordinal-1)/(1-csi)
  ui<-(m-ordinal)/(csi^2)+(ordinal-1)/((1-csi)^2)
  pri<-probcub00(m,pai,csi)
  qi<-1/(m*pri[ordinal])
  qistar<-1-(1-pai)*qi
  qitilde<-qistar*(1-qistar)
  i11<-sum((1-qi)^2)/(pai^2)
  i12<- -sum(vvi*qi*qistar)/pai
  i22<-sum(qistar*ui-(vvi^2)*qitilde)
  ####################################### Information matrix
  matinf<-matrix(c(i11,i12,i12,i22),nrow=2,byrow=T) 
  ####################################### Variance-covariance matrix
  
  
  
  if(any(is.na(matinf))==TRUE){
    warning("ATTENTION: NAs produced")
    varmat<-matrix(NA,nrow=2,ncol=2)
  } else {
    if(det(matinf)<=0){  
      warning("ATTENTION: Variance-covariance matrix NOT positive definite")
      varmat<-matrix(NA,nrow=2,ncol=2)
    } else {
      varmat<-solve(matinf)
    }
  }
  
  return(varmat)
}
