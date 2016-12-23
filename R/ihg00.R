#' @title Main function for IHG models without covariates
#' @description Estimate and validate an IHG model without covariates for given ordinal responses. 
#' @aliases ihg00
#' @usage ihg00(m, ordinal)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @return An object of the class "IHG"
#' @keywords internal 
#' @details The optimization procedure is run via "optim", option method="Brent" for constrained optimization 
#' (lower bound = 0, upper bound=1).
#' @import stats
#' @return An object of the class "IHG"

ihg00<-function(m,ordinal){
  
  tt0<-proc.time()
  freq<-tabulate(ordinal,nbins=m)
  n<-sum(freq)
  aver<-mean(ordinal)
  theta<-iniihg(m,freq) ### initial value (moment estimator of theta)
  estoptim<-optim(theta,effeihg,m=m,freq=freq,method="Brent",lower=0, upper=1,hessian=TRUE)
  theta<-estoptim$par
  errstdtheta<-1/sqrt(estoptim$hessian)
  varmat<-errstdtheta^2
  loglik<-loglikihg(m,freq,theta)
  wald2<-(theta-1/m)/errstdtheta
  pvaltheta2<-2*(1-pnorm(abs(wald2)))
  ##wald=theta/errstdtheta
  # prima: pvaltheta1=round(2*(1-pnorm(abs(wald))),20)
  # prima: pvaltheta1=round(2*(1-pnorm(abs(wald))),20)
  AICIHG<- -2*loglik+2
  BICIHG<- -2*loglik+log(n)
  csisb<-(m-aver)/(m-1); llsb<-loglikcub00(m,freq,1,csisb);
  llunif<- -n*log(m)
  nonzero<-which(freq!=0)
  logsat<- -n*log(n)+sum((freq[nonzero])*log(freq[nonzero]))
  theorpr<-probihg(m,theta)
  dissihg<-dissim(theorpr,freq/n)
  pearson<-((freq-n*theorpr))/sqrt(n*theorpr)
  X2<-sum(pearson^2)
  relares<-(freq/n-theorpr)/theorpr
  stampa<-cbind(1:m,freq/n,theorpr,pearson,relares)
  
  LL2<-1/(1+mean((freq/(n*theorpr)-1)^2))
  II2<-(loglik-llunif)/(logsat-llunif)
  FF2<-1-dissihg
  #####################################################
  durata<-proc.time()-tt0;durata<-durata[1];
  
  
  stime<-theta
  results<-list('estimates'=stime, 'loglik'=loglik,'varmat'=varmat,'BIC'=BICIHG,'time'=durata)
  
}