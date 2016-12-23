#' @title Main function for CUB models without covariates
#' @description Function to estimate and validate a CUB model without covariates for given ordinal responses.
#' @aliases cub00
#' @usage cub00(m, ordinal, maxiter, toler)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
#' @param toler Fixed error tolerance for final estimates 
#' @return An object of the class "CUB"
#' @seealso \code{\link{CUB}}, \code{\link{probbit}}, \code{\link{probcub00}}, \code{\link{loglikCUB}}
#' @keywords internal 
#' @import stats graphics
#' @examples 
#' #Running donttest option since the proposed examples require a long time run for check 
#' \donttest{
#' data(univer)
#' m=7
#' ordinal=univer[,12]
#' model=cub00(m, ordinal, maxiter=500, toler=1e-6, makeplot=TRUE,summary=TRUE)
#' parest=model$estimates  #ML estimates (pai,csi)
#' maxlik=model$loglik
#' nniter=model$niter
#' vmat=model$varmat
#' BICCUB=model$BIC
#' #############################
#' informat=univer[,8]
#' model=cub00(m,informat,maxiter=500,toler=1e-6,makeplot=FALSE,summary=FALSE)
#' parest=model$estimates   #ML estimates (pai,csi)
#' maxlik=model$loglik
#' nniter=model$niter
#' vmat=model$varmat
#' BICCUB=model$BIC
#' }

cub00<-function(m,ordinal,maxiter,toler){
  tt0<-proc.time()
  serie<-1:m
  freq<-tabulate(ordinal,nbins=m)
  n<-sum(freq)
  aver<-mean(ordinal); varcamp<-mean(ordinal^2)-aver^2;
  #######################################################
  inipaicsi<-inibest(m,freq); pai<-inipaicsi[1]; csi<-inipaicsi[2];
  ##################################################################
  loglik<-loglikcub00(m,freq,pai,csi)
  # ********************************************************************
  # ************* E-M algorithm for CUB(0,0) ***************************
  # ********************************************************************
  nniter<-1
  while(nniter<=maxiter){
    likold<-loglik
    bb<-probbit(m,csi)
    aa<-(1-pai)/(m*pai*bb)
    tau<-1/(1+aa)
    ft<-freq*tau
    averpo<-(t(serie)%*%ft)/sum(ft)
    pai<-(t(freq)%*%tau)/n  # updated pai estimate
    csi<-(m-averpo)/(m-1)   # updated csi estimate
    if(csi<0.001){
      csi<-0.001;nniter<-maxiter-1;
    }
    # print(c(pai,csi));
    loglik<-loglikcub00(m,freq,pai,csi)
    liknew<-loglik
    testll<-abs(liknew-likold) ###### print(testll); 
    # OPTIONAL printing: print(cbind(nniter,testll,pai,csi));
    if(testll<=toler) break else {loglik<-liknew}
    # OPTIONAL printing: print(loglik);
    nniter<-nniter+1
  }
  ######
  if(csi>0.999) csi<-0.99                                                             
  if(csi<0.001) csi<-0.01         ### to avoid division by 0 !!!
  if(pai<0.001) pai<-0.01         ### to avoid division by 0 !!!
  ### to ensure identifiability !!!
  ######
  AICCUB00<- -2*loglik+2*(2)
  BICCUB00<- -2*loglik+log(n)*(2)
  nomi<-rbind("pai","csi");stime<-round(c(pai,csi),5);
  ###############################
  theorpr<-probcub00(m,pai,csi)
  dissimi<-dissim(theorpr,freq/n)
  pearson<-((freq-n*theorpr))/sqrt(n*theorpr)
  X2<-sum(pearson^2)
  relares<-(freq/n-theorpr)/theorpr
  
  llunif<- -n*log(m); csisb<-(m-aver)/(m-1);
  llsb<-loglikcub00(m,freq,1,csisb)
  nonzero<-which(freq!=0)
  logsat <- -n*log(n)+sum((freq[nonzero])*log(freq[nonzero]))
  devian<-2*(logsat-loglik)
  
  LL2<-1/(1+mean((freq/(n*theorpr)-1)^2))
  II2<-(loglik-llunif)/(logsat-llunif)
  FF2<-1-dissimi
  mat1<-cbind(pai,csi,loglik,n,X2,dissimi)
  expcub<-expcub00(m,pai,csi);     varcub<-varcub00(m,pai,csi);
  
  
  varmat<-varcovcub00(m,ordinal,pai,csi)
  ### Computation of var-covar of estimates,
  if (isTRUE(varmat==matrix(NA,nrow=2,ncol=2))==TRUE){
    ddd<-matrix(NA,nrow=2,ncol=2)
    trvarmat<-ICOMP<-NA
    errstd<-wald<-pval<-rep(NA,2)
  } else {
    ddd<-diag(sqrt(1/diag(varmat)))
    nparam<-length(stime)
    trvarmat<-sum(diag(varmat))
    ICOMP<- -2*loglik + nparam*log(trvarmat/nparam) - log(det(varmat))
    errstd<-sqrt(diag(varmat));wald<-stime/errstd;
    pval<-2*(1-pnorm(abs(wald)))
    cormat<-ddd%*%varmat%*%ddd
    esq<-sqrt(varmat[1,1])/m 
  }
  ####################################################################
  # Print CUB(0,0) results of ML estimation  
#  rdiss00<-round(dissimi,3)
  stampa<-cbind(1:m,freq/n,theorpr,pearson,relares)
  durata<-proc.time()-tt0;durata<-durata[1];
  
  results<-list('estimates'=stime,'ordinal'=ordinal,'time'=durata, 
                'loglik'=loglik,'niter'=nniter,
                'varmat'=varmat,'BIC'=BICCUB00)
  
  #class(results)<-"cub"
  return(results)
}