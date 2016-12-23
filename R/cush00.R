#' @title CUSH model without covariates
#' @description Estimate and validate a CUSH model for given ordinal responses, without covariates.
#' @usage cush00(m, ordinal, shelter)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param shelter Category corresponding to the shelter choice
#' @keywords internal
#' @aliases cush00
#' @return An object of the class "GEM", "CUSH"
#' @import stats graphics


cush00<-function(m,ordinal,shelter){
  tt0<-proc.time()
  freq<-tabulate(ordinal,nbins=m); n<-length(ordinal); aver<-mean(ordinal);
  fc<-freq[shelter]/n
  deltaest<-max(0.01,(m*fc-1)/(m-1))   ### sufficient unbiased estimator
  esdelta<-sqrt((1-deltaest)*(1+(m-1)*deltaest)/(n*(m-1)))
  varmat<-esdelta^2
  wald<-deltaest/esdelta
  loglik<-loglikcush00(m,ordinal,deltaest,shelter)
  
  AICCUSH<- -2*loglik+2
  BICCUSH<- -2*loglik+log(n)
  llunif<- -n*log(m); csisb<-(m-aver)/(m-1);
  llsb<-loglikcub00(m,freq,1,csisb)
  nonzero<-which(freq!=0)
  logsat<- -n*log(n)+sum((freq[nonzero])*log(freq[nonzero]))
  devian<-2*(logsat-loglik)
  LRT<-2*(loglik-llunif)
  theorpr<-deltaest*ifelse(seq(1,m)==shelter,1,0)+(1-deltaest)/m
  pearson<-((freq-n*theorpr))/sqrt(n*theorpr)
  X2<-sum(pearson^2)
  relares<-(freq/n-theorpr)/theorpr
  diss00<-dissim(theorpr,freq/n)
  FF2<-1-diss00
  LL2<-1/(1+mean((freq/(n*theorpr)-1)^2))
  II2<-(loglik-llunif)/(logsat-llunif)  
  stampa<-cbind(1:m,freq/n,theorpr,pearson,relares)
  durata<-proc.time()-tt0; durata<-durata[1];
  ###########
  
  ####################################
  results<-list('estimates'=deltaest, 'loglik'=loglik,
                'varmat'=varmat,'BIC'= BICCUSH,'time'=durata)
}

