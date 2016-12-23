#' @title Main function for IHG models 
#' @description Main function to estimate and validate an Inverse Hypergeometric model, without or 
#' with covariates for explaining the preference parameter.
#' @aliases IHG
#' @usage IHG(Formula, data, ...)
#' @param Formula Object of class Formula.
#' @param data Data frame from which model matrices and response variables are taken.
#' @param ... Additional arguments to pass to the fitting procedure. Argument U specifies the matrix of 
#' subjects covariates to include in the model for explaining the preference parameter (not including intercept). 
#' @return An object of the class "IHG" is a list containing the following results: 
#' \item{estimates}{Maximum likelihood parameters estimates}
#' \item{loglik}{Log-likelihood function at the final estimates}
#' \item{varmat}{Variance-covariance matrix of final estimates. If no covariate is included in the model, 
#' it returns the square of the estimated standard error for the preference parameter \eqn{\theta}}
#' \item{BIC}{BIC index for the estimated model}
#' @details This is the main function for IHG models (that are nested into CUBE models, see the references below),
#'  calling for the corresponding function whenever covariates are specified. \cr
#' The parameter \eqn{\theta} represents the probability of observing a rating corresponding to the first 
#' category and is therefore a direct measure of preference, attraction, pleasantness toward the investigated item.
#'  This is reason why \eqn{\theta} is customarily referred to as the preference parameter of the IHG model.\cr
#' The estimation procedure is not iterative, so a null result for IHG$niter is produced. \cr
#' The optimization procedure is run via "optim". The variance-covariance matrix (or the estimated standard error for
#'  theta if no covariate is included) is computed as the inverse of the returned numerically differentiated 
#'  Hessian matrix (option: hessian=TRUE as argument for optim). If not positive definite,
#'  it returns a warning message and produces a matrix with NA entries.
#' @references 
#' D'Elia A. (2003). Modelling ranks using the inverse hypergeometric distribution, 
#' \emph{Statistical Modelling: an International Journal}, \bold{3}, 65--78 \cr
#' Iannario M. (2012). CUBE models for interpreting ordered categorical data with overdispersion,
#'  \emph{Quaderni di Statistica}, \bold{14}, 137--140
#' @seealso \code{\link{probihg}},  \code{\link{iniihg}}, \code{\link{loglikIHG}} 
#' @keywords internal 
#' @examples 
#' \donttest{
#' data(relgoods)
#' m<-10
#' ordinal<-na.omit(relgoods[,41]) 
#' model<-IHG(ordinal)
#' theta<-model$estimates      # ML estimates for the preference parameter theta
#' maxlik<-model$loglik
#' sqerrst<-model$varmat         # Squared standard error of theta
#' BICIHG<-model$BIC
#' #################################
#' ordinal<-relgoods[,41]
#' gender<-relgoods[,9]
#' data<-na.omit(cbind(ordinal,gender))
#' modelcov<-IHG(data[,1],U=data[,2])
#' omega<-modelcov$estimates     #  ML estimates (intercept term: omega[1])
#' maxlik<-modelcov$loglik
#' varmat<-modelcov$varmat
#' BICcov<-modelcov$BIC
#' }




IHG<-function(Formula,data,...){
  
  ellipsis.arg<-list(...)
  
  # U<-ellipsis.arg$U
  mf<-model.frame(Formula,data=data)
  ordinal<-as.numeric(model.response(mf))
  
  covtheta<-model.matrix(Formula,data=data,rhs=1)
  
  if (ncol(covtheta)==0){
    U<-NULL
  } else {
    U<-covtheta[,-1]
  }
  
  lev <- levels(factor(ordinal,ordered=TRUE))
  m <- length(lev) 
  
  if (is.null(U)){
    mod<- ihg00(m,ordinal)
  } else {
    U<-as.matrix(U)
    mod<- ihgcov(m,ordinal,U)
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




