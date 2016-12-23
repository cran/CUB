#' @title Main function for CUB models with covariates for all the components
#' @description Function to estimate and validate a CUB model for given ordinal responses, with covariates for
#'  explaining all the components and the shelter effect.
#' @aliases gecubpqs
#' @usage gecubpqs(ordinal,Y,W,X,shelter,theta0,maxiter,toler)
#' @param ordinal Vector of ordinal responses
#' @param Y Matrix of selected covariates for explaining the uncertainty component, not including intercept
#' @param W Matrix of selected covariates for explaining the feeling component, not including intercept
#' @param X Matrix of selected covariates for explaining the shelter effect, not including intercept
#' @param shelter Category corresponding to the shelter choice
#' @param theta0 Starting values for parameters explaining the shelter effect
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
#' @param toler Fixed error tolerance for final estimates 
#' @import stats graphics
#' @return An object of the class "CUB"
#' @references
#'  Piccolo D. and D'Elia A. (2008), A new approach for modelling consumers' preferences,
#'   \emph{Food Quality and Preference}, \bold{18}, 247--259 \cr
#' @references 
#' Iannario M. and Piccolo D. (2010), A new statistical model for the analysis of customer
#'  satisfaction, #' \emph{Quality Technology and Quantity management}, \bold{7}(2) 149--168 \cr
#' Iannario M. and Piccolo D. (2012), CUB models: Statistical methods and empirical evidence, in:
#' Kenett R. S. and Salini S. (eds.), \emph{Modern Analysis of Customer Surveys: with applications using R},
#'  J. Wiley and Sons, Chichester, 231--258. 
#' @keywords internal 

## EM algo for GECUB, #theta contains initial values for parameters
gecubpqs<-function(ordinal,Y,W,X,shelter,theta0,maxiter,toler){
  
  Y<-as.matrix(Y); W<-as.matrix(W);X<-as.matrix(X)
  s<-NCOL(X); q<-NCOL(W); p<-NCOL(Y);
  if (ncol(Y)==1){
    Y<-as.numeric(Y)
  }
  if (ncol(X)==1){
    X<-as.numeric(X)
  }
  if (ncol(W)==1){
    W<-as.numeric(W)
  }
  
  
  tt0<-proc.time();
  ndati<-NROW(ordinal);
  m<-length(levels(factor(ordinal,ordered=TRUE)))
  #   ****************************************************************************
  #     ************* E-M algorithm for GECUB *************************************
  #     **************************************************************************** 
  #first step
  # if (missing(theta0)){ 
  #   omega<-rep(0.1,s+1); bet<-rep(0.1,p+1); gama<-inibestgama(m,ordinal,W);
  #   #gama=rep(0.1,q+1); #In alternativa 
  #   } else  {
  bet<-theta0[1:(p+1)]; gama<-theta0[(p+2):(p+q+2)]; omega<-theta0[(p+q+3):(p+q+s+3)];
  # }
  
  loglik<-ellegecub(ordinal,Y,W,X,bet,gama,omega,shelter);
  
  psi<-c(omega,bet) #psi=omega|bet;
  
  # starting iterations
  nniter<-1;
  while(nniter <= maxiter){
    likold<-loglik;
    ### STEP 1: */
    alpha1<-logis(X,omega);        
    alpha2<-(1-alpha1)*(logis(Y,bet)); 
    alpha3<-1-alpha1-alpha2; 
    
    #/* ### STEP 2:*/
    prob1<-ifelse(as.numeric(ordinal)==shelter,1,0)   
    prob2<-bitgama(m,factor(ordinal,ordered=TRUE),W,gama)
    prob3<-1/m;  
    
    #/* ### STEP 3:*/
    num1<-alpha1*prob1;
    num2<-alpha2*prob2;
    num3<-alpha3*prob3; 
    den<-num1+num2+num3;
    
    #/* ### STEP 4:*/
    ttau1<-num1/den;                  
    ttau2<-num2/den;
    ttau3<-1-ttau1-ttau2;  
    
    
    ### STEP 5-6-7:*/
    
    datiuno<-cbind(ttau1,ttau2,X,Y)  #fissiuno=ttau1~ttau2~X~Y;
    datidue<-cbind(ttau2,ordinal,W)             #fissidue=ttau2~ordinal~W;
    
    
    param<-psi
    psioptim<-optim(param, Qunogecub, datiuno=datiuno,s=s,control=list(maxit=2000))    
    
    param<-gama;
    gamaoptim<-optim(param,Q2gecub,datidue=datidue)
    
    psinew<-psioptim$par;  
    gama<-gamaoptim$par;
    omega<-psinew[1:(s+1)]; bet<-psinew[(s+2):(p+s+2)];
    
    loglik<-ellegecub(ordinal,Y,W,X,bet,gama,omega,shelter);
    liknew<-loglik;
    testll<-abs(liknew-likold);
    if (testll <= toler) break else {likold=liknew}; 
    nniter<-nniter+1;
  }
  stime<-c(bet,gama,omega)
  #output
  loglik<-liknew;
  n<-ndati
  ####################################################################
  np<-s+p+q+3
  AICGECUB<- -2*loglik+2*np;
  BICGECUB<- -2*loglik+log(n)*np;
  
  varmat<-varcovgecub(ordinal,Y,W,X,bet,gama,omega,shelter); 
  if (isTRUE(varmat==matrix(NA,nrow=np,ncol=np))==TRUE){
    ddd<-matrix(NA,nrow=np,ncol=np)
    trvarmat<-ICOMP<-NA
    errstd<-wald<-pval<-rep(NA,np)
  } else {
    ddd<-diag(sqrt(1/diag(varmat)))
    nparam<-length(stime)
    trvarmat<-sum(diag(varmat))
    ICOMP<- -2*loglik + nparam*log(trvarmat/nparam) - log(det(varmat))
    errstd<-sqrt(diag(varmat));wald<-stime/errstd
    pval<-2*(1-pnorm(abs(wald)))
    cormat<-ddd%*%varmat%*%ddd 
  }
  
  
  durata<-proc.time()-tt0;durata<-durata[1];
  results<-list('estimates'=stime,'ordinal'=ordinal,'time'=durata, 
                'loglik'=loglik,'niter'=nniter,
                'varmat'=varmat,'BIC'=BICGECUB)
  
  
}
