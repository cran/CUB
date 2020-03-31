#' @title Main function for CUB models 
#' @description Main function to estimate and validate a CUB model for explaining uncertainty 
#' and feeling for given ratings, with or without covariates and shelter effect.
#' @aliases CUB
#' @usage CUB(Formula, data, ...)
#' @param Formula Object of class Formula.
#' @param data Data frame from which model matrices and response variables are taken.
#' @param ... Additional arguments to be passed for the specification of the model, including covariates matrices Y, W, X for 
#' #' for uncertainty, feeling and shelter, respectively.
#' @return An object of the class "GEM"-CUB": returns a list containing the following results: 
#' \item{estimates}{Maximum likelihood estimates: \eqn{(\pi, \xi)}}
#' \item{loglik}{Log-likelihood function at the final estimates}
#' \item{varmat}{Variance-covariance matrix of final estimates}
#' \item{niter}{Number of executed iterations}
#' \item{BIC}{BIC index for the estimated model}
#' @details This is the main function for CUB models, which calls for the corresponding functions whenever 
#' covariates or shelter effect are specified. It performs maximum likelihood estimation via the E-M algorithm 
#' for CUB models and extensions. The optimization procedure is run via "optim".\cr
#' It is possible to fit data with CUB models, with or without covariates 
#' for the parameters of the mixture model, and CUB models with shelter effect with no covariate included 
#' in the model. The program also checks if the estimated variance-covariance matrix is positive definite: 
#' if not, it prints a warning message and returns a matrix and related results with NA entries.
#' @references 
#' Piccolo D. and D'Elia A. (2008). A new approach for modelling consumers' preferences, \emph{Food Quality and Preference},
#' \bold{18}, 247--259 \cr
#' Iannario M. and Piccolo D. (2012). CUB models: Statistical methods and empirical evidence, in: 
#' Kenett R. S. and Salini S. (eds.), \emph{Modern Analysis of Customer Surveys: with applications using R}, 
#' J. Wiley and Sons, Chichester, 231--258\cr
#' Iannario M. (2012). Modelling \emph{shelter} choices in a class of mixture models for ordinal responses,  
#' \emph{Statistical Methods and Applications}, \bold{21}, 1--22 \cr
#' Iannario M. and Piccolo D. (2014). Inference for CUB models: a program in R, \emph{Statistica & Applicazioni}, 
#' \bold{XII} n.2, 177--204 \cr
#' Iannario M. (2016). Testing the overdispersion parameter in CUBE models,
#'  \emph{Communications in Statistics: Simulation and Computation}, \bold{45}(5), 1621--1635
#' @seealso \code{\link{probcub00}}, \code{\link{probcubp0}}, \code{\link{probcub0q}}, \code{\link{probcubpq}},
#' \code{\link{probcubshe1}}, \code{\link{loglikCUB}}, \code{\link{varmatCUB}} 
#' @keywords internal 


CUB<-function(Formula, data, ...){
  call <- match.call()
  
  ellipsis.arg<-list(...)
#   m<-ellipsis.arg[[1]]$m

#   
 # print(ellipsis.arg)
  mf<-model.frame(Formula,data=data,na.action=na.omit)
  ordinal<-as.numeric(model.response(mf))
  #m<-length(levels(factor(ordinal,ordered=TRUE)))
  
  #generico:
  #formula= ordinal~ covpai | covcsi | covshe 
  
  #solo con covariate per feeling...
  # formula=ordinal ~ 0 | covcsi| 0  
  
  covpai<-model.matrix(Formula,data=mf,rhs=1)
  covcsi<-model.matrix(Formula,data=mf,rhs=2)
  covshe<-model.matrix(Formula,data=mf,rhs=3)
  
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
  if (ncol(covshe)==0){
    X<-NULL
  } else {
    X<-covshe[,-1]
  }
  
  lista<-ellipsis.arg[[1]]

  m<-lista[['m']]
  maxiter<-lista[['maxiter']]
  toler<-lista$toler
  shelter<-lista$shelter

  
  if(!is.null(shelter)){
   
    if(m <= 4) stop("Number of ordered categories should be at least 5")
    
    if (is.null(Y) & is.null(W) & is.null(X)){
      shelter<-as.numeric(shelter)
      # shelter<-cl$shelter
      mod<-cubshe(m,ordinal,shelter,maxiter,toler)
    } else if (!is.null(Y) & !is.null(W) & !is.null(X)){
      W<-as.matrix(W)
      Y<-as.matrix(Y)
      X<-as.matrix(X)
      theta0<-lista$theta0
      s=NCOL(X); p=NCOL(Y); q=NCOL(W);
      if (is.null(theta0)){
        freq<-tabulate(ordinal,nbins=m)
        inipaicsi=inibest(m,freq)
        pai=inipaicsi[1]; bet0=log(pai/(1-pai)); bet=c(bet0,rep(0,p))
        gama<-inibestgama(m,factor(ordinal,ordered=TRUE),W)
        omega=rep(0.1,s+1);
        theta0=c(bet,gama,omega);
      }
      mod<-gecubpqs(ordinal,Y,W,X,shelter,theta0,maxiter,toler)
    } else{
      cat("CUB model with shelter effect not available for this variables specification")
    }
  } else{
    if(is.null(Y) & is.null(W)) {
      if(m <= 3) stop("Number of ordered categories should be at least 4")
      mod<-cub00(m,ordinal,maxiter,toler)
      
    }
    
    else{
      if(!is.null(Y) & is.null(W)) {
        Y<-as.matrix(Y)
        mod<-cubp0(m,ordinal,Y,maxiter,toler)
        
      }
      else{
        if(is.null(Y) & !is.null(W)) {
          W<-as.matrix(W)
          mod<-cub0q(m,ordinal,W,maxiter,toler)
          
        }
        else{
          if(!is.null(Y) & !is.null(W)) {
            Y<-as.matrix(Y)
            W<-as.matrix(W)
            mod<-cubpq(m,ordinal,Y,W,maxiter,toler)
            
          }
          else cat("Wrong variables specification")
        }
      }                            
    }
  }
  
  
  stime<-mod$estimates
  durata<-mod$time
  loglik<-as.numeric(mod$loglik)
  niter<-mod$niter
  varmat<-mod$varmat
  BIC<-as.numeric(mod$BIC)
  time<-mod$durata
  ordinal<-factor(mod$ordinal,ordered=TRUE)
  
  
  results<-list('estimates'=stime,'ordinal'=ordinal,'time'=durata,
                'loglik'=loglik,'niter'=niter,'varmat'=varmat,
                'BIC'=BIC)
  #class(results)<-"cub"
  return(results)
  
}




#################################

