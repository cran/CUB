#' @title Probability distribution of a GeCUB model
#' @aliases probgecub
#' @description Compute the probability distribution of a GeCUB model, that is a CUB model with 
#' shelter effect with covariates specified for all component.
#' @export probgecub
#' @usage probgecub(ordinal,Y,W,X,bet,gama,omega,shelter)
#' @keywords distribution
#' @param ordinal Vector of ordinal responses
#' @param Y Matrix of covariates for explaining the uncertainty component
#' @param W Matrix of covariates for explaining the feeling component
#' @param X Matrix of covariates for explaining the shelter effect
#' @param bet Vector of parameters for the uncertainty component, whose length equals 
#' NCOL(Y)+1 to include an intercept term in the model (first entry)
#' @param gama Vector of parameters for the feeling component, whose length equals 
#' NCOL(W)+1 to include an intercept term in the model (first entry)
#' @param omega Vector of parameters for the shelter effect, whose length equals 
#' NCOL(X)+1 to include an intercept term in the model (first entry)
#' @param shelter Category corresponding to the shelter choice
#' @return A vector of the same length as \code{ordinal}, whose i-th component is the
#' probability of the i-th observation according to a GeCUB model with the corresponding values 
#' of the covariates for all the components and coefficients specified in \code{bet}, \code{gama}, \code{omega}.
#' @references 
#' Iannario M. and Piccolo D. (2016b). A generalized framework for modelling ordinal data. 
#'  \emph{Statistical Methods and Applications}, \bold{25}, 163--189.\cr 

probgecub<-function(ordinal,Y,W,X,bet,gama,omega,shelter){
  
  if (is.factor(ordinal)){
    ordinal<-unclass(ordinal)
  }
  
  X<-as.matrix(X); Y<-as.matrix(Y); W<-as.matrix(W)
  
  if (ncol(W)==1){
    W<-as.numeric(W)
  }
  if (ncol(Y)==1){
    Y<-as.numeric(Y)
  }
  if (ncol(X)==1){
    X<-as.numeric(X)
  }
  
  alpha1<-logis(X,omega);
  alpha2<-(1-alpha1)*(logis(Y,bet));
  pshe<-ifelse(as.numeric(ordinal)==shelter,1,0)
  ord<-factor(ordinal,ordered=TRUE)
  m<-length(levels(ord))
  vettore<-alpha1*pshe + alpha2*(bitgama(m,ordinal,W,gama)) + (1-alpha1-alpha2)*(1/m);
  return(vettore)
}