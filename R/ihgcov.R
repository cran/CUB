#' @title Main function for IHG models with covariates
#' @description Estimate and validate an IHG model for given ordinal responses, with covariates to 
#' explain the preference parameter. 
#' @aliases ihgcov
#' @usage ihgcov(m, ordinal, U)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param U Matrix of selected covariates for the preference parameter 
#' @return An object of the class "IHG"
#' @keywords internal 
#' @details The optimization procedure is run via "optim", option method="Brent" for constrained optimization 
#' (lower bound = 0, upper bound=1).
#' @import stats graphics
#' @return An object of the class "IHG"
#' 
ihgcov<-function(m,ordinal,U){ 
  tt0<-proc.time()
  freq<-tabulate(ordinal,nbins=m); n<-length(ordinal);
  U<-as.matrix(U)
  
  if (ncol(U)==1){
    U<-as.numeric(U)
  }
  theta<-iniihg(m,freq)
  ncovar<-NCOL(U)
  nuzero<-log(theta/(1-theta))         ### initial estimate of nu_0
  nuinit<-c(nuzero,rep(0.1,ncovar)) ### initial estimate of nu vector
  ### maximize w.r.t. nu
  nu<-nuinit
  optimnu<-optim(nu,effeihgcov,ordinal=ordinal,U=U,m=m,hessian=TRUE)
  #################################################################
  # Computation of estimates and log-likelihood
  #################################################################
  nuest<-optimnu$par                #nu estimates
  loglik<-loglikihgcov(m,ordinal,U,nuest)
  HHH<-optimnu$hessian
  nparam<-length(nuest)
  
  if (det(HHH)<=0){
    warning("Variance-covariance matrix not-positive definite")
    varmat<-ddd<-matrix(NA,nrow=nparam,ncol=nparam)
    errst<-wald<-rep(NA,nparam)
    trvarmat<-ICOMP<-NA    
  } else {
    varmat<-solve(HHH)
    errst<-sqrt(diag(varmat))       ### vector
    ddd<-diag(sqrt(1/diag(varmat))) ### matrix
    wald<-nuest/errst
    trvarmat<-sum(diag(varmat))
    ICOMP<- -2*loglik + nparam*log(trvarmat/nparam) - log(det(varmat)) 
  }
  
  AICIHGCOV<- -2*loglik+2*nparam
  BICIHGCOV<- -2*loglik+nparam*log(n)
  
  nomi<-c(paste("nu",0:(nparam-1),sep="_"))
  stime<-nuest; errstd<-errst;  wald<-wald;
  pval<-2*(1-pnorm(abs(wald)))
  durata<-proc.time()-tt0;durata<-durata[1];
  
  results<-list('estimates'=stime, 'loglik'=loglik, 'varmat'=varmat,
                'BIC'=BICIHGCOV,'time'=durata)
  
}
