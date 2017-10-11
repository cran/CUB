#' @title Auxiliary function for the log-likelihood estimation of GeCUB models.
#' @aliases Q2gecub
#' @description Define the opposite one of the two scalar functions that are maximized when running the E-M algorithm
#' for GeCUB models with covariates for feeling, uncertainty and overdispersion.
#' @keywords internal 
#' @usage Q2gecub(param,datidue)
#' @param param Vector of initial estimates of parameters for the feeling component
#' @param datidue  Auxiliary matrix 


Q2gecub<-function(param,datidue){
  gama<-param; q<-NROW(gama)-1; 
  taudue<-datidue[,1];
  ordd<-datidue[,2];
  m<-length(levels(factor(ordd,ordered=TRUE)))
  W<-datidue[,3:(q+2)];
  
  return(-sum(taudue*log(bitgama(m,ordd,W,gama)))); ### 
}