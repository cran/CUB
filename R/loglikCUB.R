#' @title Log-likelihood function for CUB models
#' @aliases loglikCUB
#' @description  Compute the log-likelihood value of a CUB model fitting given data, with or without covariates to 
#' explain the feeling and uncertainty components, or for extended CUB models with shelter effect.
#' @usage loglikCUB(ordinal,m,param,Y=0,W=0,X=0,shelter=0)
#' @export loglikCUB
#' @param ordinal Vector of ordinal responses 
#' @param m Number of ordinal categories
#' @param param Vector of parameters for the specified CUB model
#' @param Y Matrix of selected covariates to explain the uncertainty component (default: no covariate is included 
#' in the model)
#' @param W Matrix of selected covariates to explain the feeling component (default: no covariate is included 
#' in the model)
#' @param X Matrix of selected covariates to explain the shelter effect (default: no covariate is included 
#' in the model)
#' @param shelter Category corresponding to the shelter choice (default: no shelter effect is included in the 
#' model)
#' @details If no covariate is included in the model, then \code{param} should be given in the form \eqn{(\pi,\xi)}. 
#' More generally, it should have the form \eqn{(\bold{\beta,\gamma)}} where, 
#' respectively, \eqn{\bold{\beta}} and \eqn{\bold{\gamma}} are the vectors of 
#' coefficients explaining the uncertainty and the feeling components, with length NCOL(Y)+1 and
#'  NCOL(W)+1 to account for an intercept term in the first entry. When shelter effect is considered, \code{param} corresponds 
#'  to the first possibile parameterization and hence should be given as \code{(pai1,pai2,csi)}. 
#'  No missing value should be present neither
#'   for \code{ordinal} nor for covariate matrices: thus, deletion or imputation procedures should be preliminarily run.
#' @seealso  \code{\link{logLik}}
#' @keywords htest
#' @examples
#' ## Log-likelihood of a CUB model with no covariate
#' m<-9; n<-300
#' pai<-0.6; csi<-0.4
#' ordinal<-simcub(n,m,pai,csi)
#' param<-c(pai,csi)
#' loglikcub<-loglikCUB(ordinal,m,param)
#' ##################################
#' ## Log-likelihood of a CUB model with covariate for uncertainty
#' \donttest{
#' data(relgoods)
#' m<-10
#' naord<-which(is.na(relgoods$Physician))
#' nacov<-which(is.na(relgoods$Gender))
#' na<-union(naord,nacov)
#' ordinal<-relgoods$Physician[-na]; Y<-relgoods$Gender[-na]
#' bbet<-c(-0.81,0.93); ccsi<-0.2
#' param<-c(bbet,ccsi)
#' loglikcubp0<-loglikCUB(ordinal,m,param,Y=Y)
#' #######################
#' ## Log-likelihood of a CUB model with covariate for feeling
#' data(relgoods)
#' m<-10
#' naord<-which(is.na(relgoods$Physician))
#' nacov<-which(is.na(relgoods$Gender))
#' na<-union(naord,nacov)
#' ordinal<-relgoods$Physician[-na]; W<-relgoods$Gender[-na]
#' pai<-0.44; gama<-c(-0.91,-0.7)
#' param<-c(pai,gama)
#' loglikcub0q<-loglikCUB(ordinal,m,param,W=W)
#' #######################
#' ## Log-likelihood of a CUB model with covariates for both parameters
#' data(relgoods)
#' m<-10
#' naord<-which(is.na(relgoods$Walking))
#' nacovpai<-which(is.na(relgoods$Gender))
#' nacovcsi<-which(is.na(relgoods$Smoking))
#' na<-union(naord,union(nacovpai,nacovcsi))
#' ordinal<-relgoods$Walking[-na]
#' Y<-relgoods$Gender[-na]; W<-relgoods$Smoking[-na]
#' bet<-c(-0.45,-0.48); gama<-c(-0.55,-0.43)
#' param<-c(bet,gama)
#' loglikcubpq<-loglikCUB(ordinal,m,param,Y=Y,W=W)
#' #######################
#' ### Log-likelihood of a CUB model with shelter effect
#' m<-7; n<-400
#' pai<-0.7; csi<-0.16; delta<-0.15
#' shelter<-5
#' ordinal<-simcubshe(n,m,pai,csi,delta,shelter)
#' pai1<- pai*(1-delta); pai2<-1-pai1-delta
#' param<-c(pai1,pai2,csi)
#' loglik<-loglikCUB(ordinal,m,param,shelter=shelter)
#' ##############
#' ### Log-likelihood of a GeCUB
#' data(univer)
#' ordinal<-univer$officeho; Y<-W<-X<-univer$gender;
#' modelgecub<-GEM(Formula(ordinal~Y|W|X),family="cub",shelter=7,maxiter=100)
#' logLik(modelgecub)
#' param<-rep(0.1,6)
#' loglik<-loglikCUB(ordinal,m=7,param=param,shelter=7,Y=Y,W=W,X=X)
#' }
#' 


loglikCUB<-function(ordinal,m,param,Y=0,W=0,X=0,shelter=0){

  if (is.factor(ordinal)){
    ordinal<-unclass(ordinal)
  }
  
  freq<-tabulate(ordinal,nbins=m)
  
  ry<-NROW(Y);   rw<-NROW(W); rx<-NROW(X); shelter<-as.numeric(shelter)
  
  if(shelter!=0){
    if (ry==1 & rw==1 & rx==1){
      pai1<-param[1]
      pai2<-param[2]
      csi<-param[3]
      loglik<-loglikcubshe(m,freq,pai1,pai2,csi,shelter)
    } else if (ry!=1 & rw !=1 & rx !=1){
      Y<-as.matrix(Y); W<-as.matrix(W);X<-as.matrix(X)
      
      if (ncol(W)==1){
        W<-as.numeric(W)
      }
      if (ncol(Y)==1){
        Y<-as.numeric(Y)
      }
      if (ncol(X)==1){
        X<-as.numeric(X)
      }
      
      
      ncy<-NCOL(Y); ncw<-NCOL(W); ncx<-NCOL(X)
      bet<-param[1:(ncy+1)]; gama<-param[(ncy+2):(ncy+ncw+2)]; omega<-param[(ncy+ncw+3):(ncy+ncw+ncx+3)];
      loglik<-ellegecub(ordinal,Y,W,X,bet,gama,omega,shelter)
    } else{
      cat("Wrong variables specification")
      loglik<-NULL
    }
    
  }else{
    
    if(ry==1 & rw==1 & rx==1) {
      pai<-param[1]
      csi<-param[2]
      loglik<-loglikcub00(m,freq,pai,csi)
    }
  
      if(ry!=1 & rw==1 & rx==1) {
        Y<-as.matrix(Y)
      
        if (ncol(Y)==1){
          Y<-as.numeric(Y)
        }
        ncy<-NCOL(Y)
       bbet<-param[1:(ncy+1)]
       ccsi<-param[length(param)]
       loglik<-loglikcubp0(m,ordinal,Y,bbet,ccsi) 
      } 
    
     if(ry==1 & rw!=1 & rx==1) {
        pai<-param[1]
        gama<-param[2:length(param)]
        W<-as.matrix(W)
        if (ncol(W)==1){
          W<-as.numeric(W)
        }
        loglik<-loglikcub0q(m,ordinal,W,pai,gama)
      }
    
     if(ry!=1 & rw!=1& rx==1) {
          ncy<-NCOL(Y)
          Y<-as.matrix(Y)
          W<-as.matrix(W)
          if (ncol(W)==1){
            W<-as.numeric(W)
          }
          if (ncol(Y)==1){
            Y<-as.numeric(Y)
          }
          bet<-param[1:(ncy+1)]
          gama<-param[(ncy+2):length(param)]
          loglik<-loglikcubpq(m,ordinal,Y,W,bet,gama)
      }
          
    }

  return(loglik)
}
