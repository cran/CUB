#' @title Log-likelihood function for CUSH models
#' @aliases loglikCUSH
#' @description  Compute the log-likelihood function for CUSH models with or without covariates 
#' to explain the shelter effect.
#' @usage loglikCUSH(ordinal,m,param,shelter,X=0)
#' @export loglikCUSH
#' @param ordinal Vector of ordinal responses (factor type)
#' @param m Number of ordinal categories
#' @param param Vector of parameters for the specified CUSH model
#' @param shelter Category corresponding to the shelter choice
#' @param X Matrix of selected covariates to explain the shelter effect (default: no covariate 
#' is included in the model)
#' @details If no covariate is included in the model, then \code{param} is the estimate of the shelter 
#' parameter (delta), otherwise \code{param} has length equal to NCOL(X) + 1 to account for an intercept  
#' term (first entry)
#' @seealso  \code{\link{GEM}}, \code{\link{logLik}}
#' @keywords htest
#' @examples
#' ## Log-likelihood of CUSH model without covariates
#' n<-300
#' m<-7
#' shelter<-2; delta<-0.4
#' ordinal<-simcush(n,m,delta,shelter)
#' loglik<-loglikCUSH(ordinal,m,param=delta,shelter)
#' #####################
#' ## Log-likelihood of CUSH model with covariates
#' data(relgoods); attach(relgoods)
#' m<-10
#' naord<-which(is.na(SocialNetwork))
#' nacov<-which(is.na(Gender))
#' na<-union(nacov,naord)
#' ordinal<-SocialNetwork[-na]; cov<-Gender[-na]
#' omega<-c(-2.29, 0.62)
#' loglikcov<-loglikCUSH(ordinal,m,param=omega,shelter=1,X=cov)

loglikCUSH<-function(ordinal,m,param,shelter,X=0){
  
  if (!is.factor(ordinal)){
    stop("Response must be an ordered factor")
  }
  
  ordinal<-unclass(ordinal)
  
  nx<-NROW(X)
  if (nx==1){
    delta<-param
    loglik<-loglikcush00(m,ordinal,delta,shelter)
    
  } else {
    omega<-param
    X<-as.matrix(X)
    
    if (ncol(X)==1){
      X<-as.numeric(X)
    }
    loglik<-loglikcushcov(m,ordinal,X,omega,shelter)  
  }
  
  return(loglik)
  
}

