#' @title CUSH model with covariates
#' @description Estimate and validate a CUSH model for ordinal responses, with covariates
#'  to explain the shelter effect.
#' @aliases cushcov
#' @usage cushcov(m, ordinal, X, shelter)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param X Matrix of selected covariates for explaining the shelter effect
#' @param shelter Category corresponding to the shelter choice
#' @keywords internal
#' @return An object of the class "GEM", "CUSH"
#' @import stats graphics


###########################################################################################
### CUSHCOV (shelter) with covariates
###########################################################################################
cushcov<-function(m,ordinal,X,shelter){ 
  tt0<-proc.time()
  freq<-tabulate(ordinal,nbins=m); n<-length(ordinal);
  fc<-freq[shelter]/n
  delta<-max(0.01,(m*fc-1)/(m-1))              ### sufficient unbiased estimator for a CUSH model
  X<-as.matrix(X)
  if (ncol(X)==1){
    X<-as.numeric(X)
  }
  
  ncovar<-NCOL(X)
  omzero<-log(delta/(1-delta))         ### initial estimate of omega_0
  omegainit<-c(omzero,rep(0.1,ncovar)) ### initial estimate of omega vector
  ### maximize w.r.t. omega 
  XX<-cbind(1,X)
  esternocush<-cbind(ordinal,XX)
  paravec<-omegainit 
  shelter<-shelter
  optimomega<-optim(paravec,effecush,esternocush,shelter=shelter,m=m,gr=NULL,hessian=TRUE)
  #################################################################
  # Computation of estimates and log-likelihood
  #################################################################
  omegaest<-optimomega$par               #omega estimates
  loglik<-loglikcushcov(m,ordinal,X,omegaest,shelter) #loglik at the maximum
  HHH<-optimomega$hessian
  nparam<-length(omegaest)
  
  if (det(HHH)<=0){
    warning("Variance-Covariance matrix is not positive definite")
    varmat<-ddd<-cormat<-matrix(NA,nrow=nparam,ncol=nparam)
    trvarmat<-ICOMP<-NA
    errst<-wald<-pval<-rep(NA,nparam)  
  } else {
    varmat<-solve(HHH)
    errst<-sqrt(diag(varmat))       ### vector
    ddd<-diag(sqrt(1/diag(varmat))) ### matrix
    wald<-omegaest/errst
    trvarmat<-sum(diag(varmat))
    ICOMP<- -2*loglik + nparam*log(trvarmat/nparam) - log(det(varmat)) ## added
  }

  AICCUSH<- -2*loglik+2*nparam
  BICCUSH<- -2*loglik+nparam*log(n)
  
  nomi<-c(paste("omega",0:(nparam-1),sep="_"))
  stime<-omegaest; errstd<-errst;  wald<-wald;
  pval<-2*(1-pnorm(abs(wald)))
  durata<-proc.time()-tt0;durata<-durata[1];

  results<-list('estimates'=stime, 'loglik'=loglik, 'varmat'=varmat,'BIC'=BICCUSH,
                'time'=durata)
  
}

