#' @title Main function for CUB models with a shelter effect
#' @description Estimate and validate a CUB model with a shelter effect.
#' @aliases cubshe
#' @usage cubshe(m, ordinal, shelter, maxiter, toler)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param shelter Category corresponding to the shelter choice
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
#' @param toler Fixed error tolerance for final estimates 
#' @param makeplot Logical: if TRUE (default), the algorithm returns a graphical plot comparing fitted 
#' probabilities and observed relative frequencies
#' @param summary Logical: if TRUE, summary results of the fitting procedure are displayed on screen
#' @return An object of the class "CUB"
#' @import stats
#' @references
#' Iannario M. (2012). Modelling \emph{shelter} choices in a class of mixture models for ordinal responses,  
#' \emph{Statistical Methods and Applications}, \bold{21}, 1--22
#' @keywords internal #models


cubshe<-function(m,ordinal,shelter,maxiter,toler){
  tt0<-proc.time()
  
  #######
  ###########
  
  serie<-1:m; freq<-tabulate(ordinal,nbins=m); n<-sum(freq); 
  
  ##########
  ###########
  
  aver<-mean(ordinal); varcamp<-mean(ordinal^2)-aver^2;
  dd<-ifelse(serie==shelter,1,0)
  #####################################################################
  vett<-inibest(m,freq)
  pai1<-vett[1]; csi<-vett[2]
  

  fc<-freq[shelter]/n
  
  deltaini<-max(0.01,(m*fc-1)/(m-1)) #deltaini=runif(1,0,0.3)  
  pai2<-max(0.01,1-deltaini-pai1)
  ################################################################
  loglik<-loglikcubshe(m,freq,pai1,pai2,csi,shelter)
  

  # ********************************************************************
  # ************* E-M algorithm for CUBSHE *****************************
  # ********************************************************************
  nniter<-1
  while(nniter<=maxiter){
    
    #############

    likold<-loglik
    bb<-probbit(m,csi)
    tau1<-pai1*bb
    tau2<-pai2*(1/m)
    denom<-tau1+tau2+(1-pai1-pai2)*dd
    tau1<-tau1/denom
    tau2<-tau2/denom        
    tau3<-1-tau1-tau2
    numaver<-sum(serie*freq*tau1)
    denaver<-sum(freq*tau1)
    averpo<-numaver/denaver
    pai1<-sum(freq*tau1)/n    #updated pai1 estimate
    pai2<-sum(freq*tau2)/n    #updated pai2 estimate
    csi<-(m-averpo)/(m-1)     #updated csi estimate
    
    
    if(csi<0.001){
      csi<-0.001;nniter<-maxiter-1;
    }
    
    loglik<-loglikcubshe(m,freq,pai1,pai2,csi,shelter)
    liknew<-loglik
    testll<-abs(liknew-likold)
    #print(cbind(nniter,testll,pai1,pai2,csi,loglik));
    if(testll<=toler) break else {loglik<-liknew}
    nniter<-nniter+1
  }
  
  ######
  if(csi>0.999) csi<-0.99   ###????????### to avoid division by 0 !!!
  if(csi<0.001) csi<-0.01  ###????### to avoid division by 0 !!!
  if(pai1<0.001) pai1<-0.01   ###?????         ### to ensure identifiability !!!
  
  ####################################################################
  # Compute asymptotic standard errors of ML estimates
  ####################################################################
  varmat<-varcovcubshe(m,pai1,pai2,csi,shelter,n)
  nomi<-rbind("pai1","pai2","csi")
  stime<-c(pai1,pai2,csi)
  nparam<-length(stime)
  delta<-1-pai1-pai2
  paistar<-pai1/(pai1+pai2)
  if (isTRUE(varmat==matrix(NA,nrow=nparam,ncol=nparam))==TRUE){
    ddd<-cormat<-matrix(NA,nrow=nparam,ncol=nparam)
    ICOMP<-trvarmat<-NA
    esdelta<-pvaldelta<-espaistar<-pvalpaistar<-NA
    errstd<-wald<-pval<-rep(NA,nparam)
  } else {
    ddd<-diag(sqrt(1/diag(varmat)))
    esdelta<-sqrt(varmat[1,1]+varmat[2,2]+2*varmat[1,2])
    pvaldelta<-round(2*(1-pnorm(abs(delta/esdelta))),20)
    espaistar<-sqrt((pai1^2*varmat[2,2]+pai2^2*varmat[1,1]-2*pai1*pai2*varmat[1,2]))/(pai1+pai2)^2
    pvalpaistar<-round(2*(1-pnorm(abs(paistar/espaistar))),20)  
    trvarmat<-sum(diag(varmat))
    ICOMP<- -2*loglik + nparam*log(trvarmat/nparam) - log(det(varmat))
    errstd<-sqrt(diag(varmat))  
    wald<-stime/errstd
    pval<-round(2*(1-pnorm(abs(wald))),20)
  }
  theorpr<-probcubshe1(m,pai1,pai2,csi,shelter)
  dissshe<-dissim(theorpr,freq/n)
  llunif<- -n*log(m); csisb<-(m-aver)/(m-1);
  llsb<-loglikcub00(m,freq,1,csisb)
  llunif<- -n*log(m)
  nonzero<-which(freq!=0)
  logsat<- -n*log(n)+sum((freq[nonzero])*log(freq[nonzero]))
  
  pearson<-((freq-n*theorpr))/sqrt(n*theorpr)
  X2<-sum(pearson^2)
  relares<-(freq/n-theorpr)/theorpr
  stampa<-cbind(1:m,freq/n,theorpr,pearson,relares)
  
  LL2<-1/(1+mean((freq/(n*theorpr)-1)^2))
  II2<-(loglik-llunif)/(logsat-llunif)
  FF2<-1-dissshe
  AICCUBshe<- -2*loglik+2*(3)
  BICCUBshe<- -2*loglik+log(n)*(3)
  durata<-proc.time()-tt0;durata<-durata[1];
  # 
  results<-list('estimates'=stime,'ordinal'=ordinal,'time'=durata,
                'loglik'=loglik,'niter'=nniter,'varmat'=varmat,
                'BIC'=BICCUBshe)
  #class(results)<-"cub"
  return(results)
}