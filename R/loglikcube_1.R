#' @title Log-likelihood function for CUBE models
#' @aliases loglikCUBE
#' @description  Compute the log-likelihood function for CUBE models. It is possible to include 
#'  covariates in the model for explaining the feeling component or all the three parameters.
#' @usage loglikCUBE(ordinal,m,param,Y=0,W=0,Z=0)
#' @export loglikCUBE
#' @param ordinal Vector of ordinal responses (factor type)
#' @param m Number of ordinal categories
#' @param param Vector of parameters for the specified CUBE model
#' @param Y Matrix of selected covariates to explain the uncertainty component (default: no covariate is included 
#' in the model)
#' @param W Matrix of selected covariates to explain the feeling component (default: no covariate is included 
#' in the model)
#' @param Z Matrix of selected covariates to explain the overdispersion component (default: no covariate is included 
#' in the model)
#' @details If no covariate is included in the model, then \code{param} has the form \eqn{(\pi,\xi,\phi)}. More generally, 
#' it has the form \eqn{(\bold{\beta,\gamma,\alpha)}} where, respectively, \eqn{\bold{\beta}},\eqn{\bold{\gamma}}, \eqn{\bold{\alpha}}
#'  are the vectors of  coefficients explaining the uncertainty, the feeling and the overdispersion components, with length NCOL(Y)+1, 
#'  NCOL(W)+1, NCOL(Z)+1 to account for an intercept term in the first entry.
#' @seealso  \code{\link{logLik}}
#' @keywords htest
#' @examples
#' #### Log-likelihood of a CUBE model with no covariate
#' m<-7; n<-400
#' pai<-0.83; csi<-0.19; phi<-0.045
#' ordinal<-simcube(n,m,pai,csi,phi)
#' loglik<-loglikCUBE(ordinal,m,param=c(pai,csi,phi))
#' ##################################
#' #### Log-likelihood of a CUBE model with covariate for feeling
#' data(relgoods)
#' attach(relgoods)
#' m<-10
#' nacov<-which(is.na(BirthYear))
#' naord<-which(is.na(Tv))
#' na<-union(nacov,naord)
#' age<-2014-BirthYear[-na]
#' lage<-log(age)-mean(log(age))
#' ordinal<-Tv[-na]; W<-lage
#' pai<-0.63; gama<-c(-0.61,-0.31); phi<-0.16
#' param<-c(pai,gama,phi)
#' loglik<-loglikCUBE(ordinal,m,param,W=W)
#' ########## Log-likelihood of a CUBE model with covariates for all parameters
#' Y<-W<-Z<-lage
#' bet<-c(0.18, 1.03); gama<-c(-0.6, -0.3); alpha<-c(-2.3,0.92)
#' param<-c(bet,gama,alpha)
#' loglik<-loglikCUBE(ordinal,m,param,Y=Y,W=W,Z=Z)


loglikCUBE <-
function(ordinal,m,param,Y=0,W=0,Z=0){
  
  if (!is.factor(ordinal)){
    stop("Response must be an ordered factor")
  }
  
  ordinal<-unclass(ordinal)
  
  ry<-NROW(Y); rw<-NROW(W); rz<-NROW(Z);
  
  freq<-tabulate(ordinal,nbins=m)
  if(ry==1 & rw==1 & rz==1) {
    pai<-param[1]; csi<-param[2]; phi<-param[3];
    loglik<-loglikcube(m,freq,pai,csi,phi) 
  } else if (ry==1 & rz==1 & rw >1){
     ncw<-NCOL(W)
     W<-as.matrix(W)
     if (ncol(W)==1){
       W<-as.numeric(W)
     }
     pai<-param[1]; gama<-param[2:(ncw+2)]; phi<-param[length(param)];
     loglik<-loglikcubecsi(m,ordinal,W,pai,gama,phi)
  } else if(ry>1 & rz>1 & rw >1){
    Y<-as.matrix(Y);  W<-as.matrix(W);    Z<-as.matrix(Z)
    if (ncol(W)==1){
      W<-as.numeric(W)
    }
    if (ncol(Y)==1){
      Y<-as.numeric(Y)
    }
    if (ncol(Z)==1){
      Z<-as.numeric(Z)
    }
   ncy<-NCOL(Y)
   ncw<-NCOL(W)
   bet<-param[1:(ncy+1)]; gama<-param[(ncy+2):(ncy+ncw+2)]; alpha<-param[(ncy+ncw+3):length(param)];
    loglik<-loglikcubecov(m,ordinal,Y,W,Z,bet,gama,alpha)
 } else {
    cat("CUBE models not available for this variables specification")
  }
return(loglik)
}
