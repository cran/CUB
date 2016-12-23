#' @title Auxiliary function for the log-likelihood estimation of GeCUB models.
#' @aliases Qunogecub
#' @description Define the opposite one of the two scalar functions that are maximized when running the E-M algorithm
#' for GeCUB models with covariates for feeling, uncertainty and overdispersion.
#' @keywords internal 
#' @usage Qunogecub(param,datiuno,s)
#' @param param Vector of initial estimates of parameters for the uncertainty component
#' @param datiuno  Auxiliary matrix 
#' @param s Number of covariates to explain the shelter effect


Qunogecub<-function(param,datiuno,s){
  p<-NROW(param)-s-2; 
  omega<-param[1:(s+1)];
  bet<-param[(s+2):(p+s+2)];
  tauno<-datiuno[,1];
  taudue<-datiuno[,2];
  covar<-datiuno[,3:(s+p+2)];
  X<-covar[,1:s];        Y=covar[,(s+1):(s+p)];
  alpha1<-logis(X,omega);
  alpha2<-(1-alpha1)*logis(Y,bet); 
  esse1<-sum(tauno*log(alpha1));
  esse2<- sum(taudue*log(alpha2)); 
  esse3<- sum((1-tauno-taudue)*log(1-alpha1-alpha2))   
  return(-esse1-esse2-esse3);
}
