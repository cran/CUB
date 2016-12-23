#' @title Main function for CUSH models 
#' @description Main function to estimate and validate a CUSH model for ordinal responses, with or without covariates
#'  to explain the shelter effect.
#' @aliases CUSH
#' @usage CUSH(Formula,data,...)
#' @param Formula Object of class Formula.
#' @param data Data frame from which model matrices and response variables are taken.
#' @param ... Additional arguments to pass to the fitting procedure. Argument X specifies the matrix of 
#' subjects covariates to include in the model for explaining the shelter effect (not including intercept). 
#' @return An object of the class "CUSH" is a list containing the following results: 
#' \item{estimates}{Maximum likelihood parameters estimates}
#' \item{loglik}{Log-likelihood function at the final estimates}
#' \item{varmat}{Variance-covariance matrix of final estimates (if X=0, it returns the square of the estimated standard error 
#' for the shelter parameter \eqn{\delta})}
#' \item{BIC}{BIC index for the estimated model}
#' @details The estimation procedure is not iterative, so a null result for CUSH$niter is produced.
#' The optimization procedure is run via "optim". If covariates are included, the variance-covariance matrix 
#' is computed as the inverse of the returned numerically differentiated Hessian matrix (option: hessian=TRUE
#'  as argument for "optim"). If not positive definite, it returns a warning message and produces a matrix 
#'  with NA entries.
#' @references 
#' Capecchi S. and Piccolo D. (2015). Dealing with heterogeneity/uncertainty in sample survey with ordinal data, 
#' \emph{IFCS Proceedings, University of Bologna} \cr
#' Capecchi S. and Iannario M. (2016). Gini heterogeneity index for detecting uncertainty in ordinal data surveys,
#'  \emph{Metron} - DOI: 10.1007/s40300-016-0088-5
#' @seealso \code{\link{loglikCUSH}}
#' @keywords internal 


CUSH<-function(Formula,data,...){  
  
  ellipsis.arg<-list(...)
  
  mf<-model.frame(Formula,data=data)
  ordinal<-as.numeric(model.response(mf))
  
  #covpai<-model.matrix(Formula,data=data,rhs=1)
  #covcsi<-model.matrix(Formula,data=data,rhs=2)
  covshe<-model.matrix(Formula,data=data,rhs=1)
  
  # if (ncol(covpai)==0){
  #   Y<-NULL
  # } else {
  #   stop("Error: no uncertainty parameter modelled with CUSH model")
  #   #Y<-as.numeric(covpai[,-1])
  # }
  # if (ncol(covcsi)==0){
  #   W<-NULL
  # } else {
  #   stop("Error: no feeling parameter modelled with CUSH model")
  #   #W<-as.numeric(covcsi[,-1])
  # }
  if (ncol(covshe)==0){
    X<-NULL
  } else {
    X<-covshe[,-1]
  }
  
  lev <- levels(factor(ordinal,ordered=TRUE))
  m <- length(lev) 
  
  shelter<-ellipsis.arg$shelter
  
  if (is.null(shelter)) stop("Shelter category missing")
  
  
  if(is.null(X)){ 
    mod<- cush00(m,ordinal,shelter)
  } else {
    X<-as.matrix(X)
    mod<-cushcov(m,ordinal,X,shelter)
  }
  stime<-mod$estimates
  durata<-mod$time
  loglik<-as.numeric(mod$loglik)
  niter<-1
  varmat<-mod$varmat
  BIC<-as.numeric(mod$BIC)
  ordinal<-factor(ordinal,ordered=TRUE)
  
  results<-list('estimates'=stime,'ordinal'=ordinal,'time'=durata,
                'loglik'=loglik,'niter'=niter,'varmat'=varmat,
                'BIC'=BIC)
  return(results)
  
}