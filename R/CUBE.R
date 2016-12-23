#' @title Main function for CUBE models
#' @description Main function to estimate and validate a CUBE model for given ratings, 
#' explaining uncertainty, feeling and overdispersion.
#' @aliases CUBE
#' @usage CUBE(Formula,data,...)
#' @param Formula Object of class Formula.
#' @param data Data frame from which model matrices and response variables are taken.
#' @param ... Additional arguments to be passed for the specification of the model, Including Y, W, Z for 
#'  explanatory variables for uncertainty, feeling and overdispersion. Set expinform=TRUE if inference should
#'   be based on expected information matrix for model with no covariate. Set starting = ... to pass initial 
#'   values for EM iterations.
#' @return An object of the class "GEM"-"CUBE" is a list containing the following results: 
#' \item{estimates}{Maximum likelihood estimates: \eqn{(\pi, \xi, \phi)}}
#' \item{loglik}{Log-likelihood function at the final estimates}
#' \item{varmat}{Variance-covariance matrix of final estimates}
#' \item{niter}{Number of executed iterations}
#' \item{BIC}{BIC index for the estimated model}
#' @details It is the main function for CUBE models, calling for the corresponding functions whenever
#'  covariates are specified: it is possible to select covariates for explaining all the three parameters
#'   or only the feeling component. \cr
#' The program also checks if the estimated variance-covariance matrix is positive definite: if not,
#'  it prints a warning message and returns a matrix and related results with NA entries.
#'  The optimization procedure is run via "optim". If covariates are included only for feeling,
#' the variance-covariance matrix is computed as the inverse of the returned numerically differentiated
#'  Hessian matrix (option: hessian=TRUE as argument for "optim"), and the estimation procedure is not
#'  iterative, so a NULL result for $niter is produced.
#'  If the estimated variance-covariance matrix is not positive definite, the function returns a 
#'  warning message and produces a matrix with NA entries.
#' @references 
#' Iannario M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data, 
#' \emph{Communications in Statistics - Theory and Methods}, \bold{43}, 771--786 \cr
#' Piccolo D. (2015). Inferential issues for CUBE models with covariates,
#'  \emph{Communications in Statistics. Theory and Methods}, \bold{44}(23), 771--786. \cr
#' Iannario M. (2015). Detecting latent components in ordinal data with overdispersion by means
#'  of a mixture distribution, \emph{Quality & Quantity}, \bold{49}, 977--987 \cr
#'  Iannario M. (2016). Testing the overdispersion parameter in CUBE models. 
#'  \emph{Communications in Statistics: Simulation and Computation}, \bold{45}(5), 1621--1635.\cr
#' @seealso \code{\link{probcube}}, \code{\link{loglikCUBE}}, \code{\link{loglikcuben}},  \code{\link{inibestcube}},
#'  \code{\link{inibestcubecsi}}, \code{\link{inibestcubecov}},
#' \code{\link{varmatCUBE}}
#' @keywords internal #models
#' @examples 
#' \donttest{
#' data(relgoods)
#' ordinal<-na.omit(relgoods[,37])  
#' model<-CUBE(ordinal,starting=c(0.1,0.1,0.1))  
#' model$estimates        # Final ML estimates
#' model$loglik           # Maximum value of the log-likelihood function
#' model$varmat         
#' model$niter
#' model$BIC
#' ######################## 
#' ordinal<-relgoods[,40]
#' cov<-relgoods[,2]
#' nona<-na.omit(cbind(ordinal,cov))
#' modelcovcsi<-CUBE(nona[,1],W=nona[,2])
#' modelcov<-CUBE(nona[,1],Y=nona[,2],W=nona[,2], Z=nona[,2])
#' modelcov$BIC
#' modelcovcsi$BIC
#' #######################################
#' data(univer)
#' ordinal<-univer[,8]
#' starting<-inibestcube(m,ordinal)
#' model<-CUBE(ordinal,starting=starting)
#' }
#' 
#' 
CUBE<-function(Formula,data,...){
  
  ellipsis.arg<-list(...)
  
  mf<-model.frame(Formula,data=data)
  
  ordinal<-as.numeric(model.response(mf))
  
  #formula given as Formula=ordinal~covpai|covcsi|covphi
  covpai<-model.matrix(Formula,data=data,rhs=1)
  covcsi<-model.matrix(Formula,data=data,rhs=2)
  covphi<-model.matrix(Formula,data=data,rhs=3)
  
  if (ncol(covpai)==0){
    Y<-NULL
  } else {
    Y<-covpai[,-1]
  }
  if (ncol(covcsi)==0){
    W<-NULL
  } else {
    W<-covcsi[,-1]
  }
  if (ncol(covphi)==0){
    Z<-NULL
  } else {
    Z<-covphi[,-1]
  }
  
  lista<-ellipsis.arg[[1]]
  # Y<-lista$Y
  # W<-lista$W
  # Z<-lista$Z
  #m<-lista$m
  
  maxiter<-lista$maxiter
  toler<-lista$toler
  starting<-lista$starting
  expinform<-lista$expinform
  
  lev <- levels(factor(ordinal,ordered=TRUE))
  m <- length(lev) 
  
  # if (is.null(maxiter)){
  #  maxiter <- 1000
  # }
  # if (is.null(toler)){
  #   toler <- 1e-6
  # }
  if (is.null(expinform)){
    expinform <-  FALSE
  }
  
  if(is.null(starting)){
    if(is.null(W) & is.null(Y) & is.null(Z)) {
      starting<-inibestcube(m,factor(ordinal,ordered=TRUE))
    }else if(is.null(Y) & is.null(Z) & !is.null(W)){
      W<-as.matrix(W)
      initial<-inibestcube(m,factor(ordinal,ordered=TRUE))
      starting<-inibestcubecsi(m,factor(ordinal,ordered=TRUE),W,initial,maxiter=500,toler=1e-6)
    } else {
      W<-as.matrix(W); Y<-as.matrix(Y); Z<-as.matrix(Z);
      starting<-inibestcubecov(m,factor(ordinal,ordered=TRUE),Y,W,Z)
    }
  }
  
  
  if(is.null(W) & is.null(Y) & is.null(Z)) {
    mod<-cube000(m,ordinal,starting,maxiter,toler,expinform)
  } else if (is.null(Y) & is.null(Z) & !is.null(W)){
    W<-as.matrix(W)
    mod<-cubecsi(m,ordinal,W,starting,maxiter,toler)  
  } else if(!is.null(Y) & !is.null(W) & !is.null(Z)){
    W<-as.matrix(W); Y<-as.matrix(Y); Z<-as.matrix(Z);
    mod<-cubecov(m,ordinal,Y,W,Z,starting,maxiter,toler=1e-2)
  } else {
    cat("CUBE models not available for this variables specification")
  }
  
  
  stime<-mod$estimates
  durata<-mod$time
  loglik<-as.numeric(mod$loglik)
  niter<-mod$niter
  varmat<-mod$varmat
  BIC<-as.numeric(mod$BIC)
  ordinal<-factor(mod$ordinal,ordered=TRUE)
  
  results<-list('estimates'=stime,'ordinal'=ordinal,'time'=durata,
                'loglik'=loglik,'niter'=niter,'varmat'=varmat,
                'BIC'=BIC)
  #class(results)<-"cube"
  return(results)
  
}