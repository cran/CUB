#' @title Log-likelihood function for IHG models
#' @aliases loglikIHG
#' @description  Compute the log-likelihood function for IHG models with or without covariates 
#' to explain the preference parameter.
#' @usage loglikIHG(ordinal,m,param,U=0)
#' @export loglikIHG
#' @param ordinal Vector of ordinal responses (factor type)
#' @param m Number of ordinal categories
#' @param param Vector of parameters for the specified IHG model
#' @param U Matrix of selected covariates to explain the preference parameter (default: no covariate is included 
#' in the model)
#' @details If no covariate is included in the model, then \code{param} is the estimate of the preference
#' parameter (theta), otherwise \code{param} has length equal to NCOL(U) + 1 to account for an intercept  
#' term (first entry)
#' @seealso  \code{\link{GEM}}, \code{\link{logLik}}
#' @keywords htest
#' @examples
#' #### Log-likelihood of an IHG model with no covariate
#' m<-10; theta<-0.14; n<-300
#' ordinal<-simihg(n,m,theta)
#' loglik<-loglikIHG(ordinal,m,param=theta)
#' ##################################
#' #### Log-likelihood of a IHG model with covariate 
#' data(relgoods)
#' attach(relgoods)
#' m<-10
#' naord<-which(is.na(HandWork))
#' nacov<-which(is.na(Gender))
#' na<-union(naord,nacov)
#' ordinal<-HandWork[-na]; U<-Gender[-na]
#' nu<-c(-1.55,-0.11)     # first entry: intercept term
#' loglik<-loglikIHG(ordinal,m,param=nu,U=U)




loglikIHG <-
function(ordinal,m,param,U=0){
  if (!is.factor(ordinal)){
    stop("Response must be an ordered factor")
  }
  
  ordinal<-unclass(ordinal)

  
  nu<-NROW(U)
  if (nu==1){
    theta<-param
    freq<-tabulate(ordinal,nbins=m)
    loglik<-loglikihg(m,freq,theta)

  } else {
   nu<-param
   U<-as.matrix(U)
   if (ncol(U)==1){
     U<-as.numeric(U)
   }
 
   loglik<-loglikihgcov(m,ordinal,U,nu)  
  }

return(loglik)
}
