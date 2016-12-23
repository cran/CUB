#' @title Variance-covariance matrix for CUB models with shelter effect
#' @description Compute the variance-covariance matrix of parameter estimates of a CUB model with shelter effect.
#' @aliases varcovcubshe
#' @usage varcovcubshe(m, pai1, pai2, csi, shelter, n)
#' @param m Number of ordinal categories
#' @param pai1 Parameter of the mixture distribution: mixing coefficient for the shifted Binomial component
#' @param pai2 Second parameter of the mixture distribution: mixing coefficient for the discrete Uniform component
#' @param csi Feeling parameter
#' @param shelter Category corresponding to the shelter choice
#' @param n Number of observations
#' @seealso \code{\link{probcubshe1}}
#' @keywords internal
#' @details The function checks if the variance-covariance matrix is positive-definite: if not, it returns a warning
#'  message and produces a matrix with NA entries.
#' @references Iannario, M. (2012), Modelling shelter choices in ordinal data surveys. 
#' Statistical Modelling and Applications, \bold{21}, 1--22

varcovcubshe <-function(m,pai1,pai2,csi,shelter,n){
  pr<-probcubshe1(m,pai1,pai2,csi,shelter)
  dd<-rep(0,m);dd[shelter]<-1;
  bb<-probbit(m,csi)
  ########################
  aaa<-bb-dd
  bbb<-(1/m)-dd
  c4<-pai1*bb*(m-(1:m)-csi*(m-1))/(csi*(1-csi))
  atilde<-aaa/pr;  btilde<-bbb/pr;  ctilde<-c4/pr;
  
  d11<-sum(aaa*atilde);  d22<-sum(bbb*btilde);   dxx<-sum(c4*ctilde);
  d12<-sum(bbb*atilde);  d1x<-sum(c4*atilde);    d2x<-sum(c4*btilde);
  
  ### Information matrix 
  matinf<-matrix(c(d11,d12,d1x,d12,d22,d2x,d1x,d2x,dxx),nrow=3,byrow=T) 
  ### Var-covar matrix 
  
  
  if(any(is.na(matinf))==TRUE){
    warning("ATTENTION: NAs produced")
    varmat<-matrix(NA,nrow=3,ncol=3)
  } else {
    if(det(matinf)<=0){  
      warning("ATTENTION: Variance-covariance matrix NOT positive definite")
      varmat<-matrix(NA,nrow=3,ncol=3)
    } else {
      varmat<-solve(matinf)/n
    }
  }
  
  
  return(varmat)
}
