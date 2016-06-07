#' @title Main function for CUB models without covariates
#' @description Function to estimate and validate a CUB model without covariates for given ordinal responses.
#' @aliases cub00
#' @usage cub00(m, ordinal, maxiter, toler, makeplot,summary)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
#' @param toler Fixed error tolerance for final estimates 
#' @param makeplot Logical: if TRUE, the function returns a graphical plot comparing fitted probabilities 
#' and observed relative frequencies
#' @param summary Logical: if TRUE, summary results of the fitting procedure are displayed on screen
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


cub00<-function(m,ordinal,maxiter,toler,makeplot,summary){
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
    errstd<-round(sqrt(diag(varmat)),5);wald<-round(stime/errstd,5);
    pval<-round(2*(1-pnorm(abs(wald))),20)
    cormat<-round(ddd%*%varmat%*%ddd,5)  
    esq<-sqrt(varmat[1,1])/m 
  }
  ####################################################################
  # Print CUB(0,0) results of ML estimation  
  rdiss00<-round(dissimi,3)
  stampa<-cbind(1:m,freq/n,theorpr,pearson,relares)
  durata<-proc.time()-tt0;durata<-durata[1];
  
  ####################################################################
  if (summary==TRUE){
    cat("\n")
    cat("=======================================================================","\n")
    #cat("================== CUB Package: version 0.1 (April 2016) ==============","\n")
    #cat("=======================================================================","\n")
    cat("=====>>> C U B (0,0) model  <<<=====   ML-estimates via E-M algorithm  ","\n")
    cat("=======================================================================","\n")
    cat(" m=", m," Sample size: n=", n," Iterations=",nniter," Maxiter=",maxiter,"\n")
    cat("=======================================================================","\n")
    cat("parameters  ML-estimates  stand.errors    estimates/stand.errors       ","\n")
    cat("=======================================================================","\n")
    for(i in 1:2){
      cat(nomi[i],"        ",stime[i],"        ",errstd[i],"  ",wald[i],"      ","\n")
    }
    cat("=======================================================================","\n")
    cat("Parameters correlation matrix","\n") 
    print(round(ddd%*%varmat%*%ddd,5))
    ##############################################################################
    cat("=======================================================================","\n")
    cat("Log-lik(pai^,csi^) =",round(loglik,digits=8),"\n")
    cat("Mean Log-likelihood=",round(loglik/n,digits=8),"\n")
    cat("Log-lik(saturated) =",round(logsat,digits=8),"\n")
    cat("Deviance           =",round(devian,digits=8),"\n")
    cat("-----------------------------------------------------------------------","\n")
    cat("Log-lik(UNIFORM)         =",round(llunif,digits=8),"\n")
    cat("Log-lik(Shifted-BINOMIAL)=",round(llsb,digits=8),"\n")
    cat("-----------------------------------------------------------------------","\n")
    cat("AIC-CUB00          =",round(AICCUB00,digits=8),"\n")
    cat("BIC-CUB00          =",round(BICCUB00,digits=8),"\n")
    cat("ICOMP-CUB          =",round(ICOMP,digits=8),"\n")
    cat("=======================================================================","\n")
    cat("Pearson Fitting measure    ==>  X^2 =",X2,"(p-val.=",1-pchisq(X2,m-3),")","\n")
    cat("Lik-based fitting measure  ==>  L^2 =",LL2,"\n")
    cat("Relative Log-lik index     ==>  I   =",round(II2,digits=5),"\n")
    cat("F^2 fitting measure        ==>  F^2 =",round(FF2,digits=5),"\n")
    cat("Normed Dissimilarity index ==>  Diss=",dissimi,"\n")
    cat("=======================================================================","\n")
    cat("Observed average          =",   aver," Sample variance        =",varcamp,"\n")
    cat("Expectation of R~CUB(0,0) =", expcub," Variance of R~CUB(0,0) =",varcub,"\n")
    cat("=======================================================================","\n")
    cat("(R=r) Observed CUB-prob","Pearson","Relative res.","\n")
    print(stampa,digits=5)
    cat("=======================================================================","\n")
    cat("Elapsed time=",durata,"seconds","=====>>>",date(),"\n")
    cat("=======================================================================","\n")  
    
  }
  
  if(makeplot==TRUE){
    stringtitle<-"CUB model";
    plot(cbind(1:m,1:m),cbind(theorpr,(freq/n)),las=1,
         main=paste(stringtitle,  "     (Diss =",round(dissimi,digits=4),")"),
         xlim=c(1,m),ylim=c(0.0,1.1*max(theorpr,(freq/n))),
         xlab="Ordinal values of R=1,2,...,m",
         ylab="Observed freq.(dots) and fitted prob. (circles)");
    ###
    points(1:m,theorpr,pch=21,cex=1.5,lwd=2.0,type="b",lty=3); ### ex pch=8,col="red"
    points(1:m,freq/n,pch=16,cex=1.25,lwd=1.5);
    abline(h=0);
  }
  #   setClass("cub", representation(estimates = "numeric", loglik = "numeric",niter="numeric",varmat="matrix",BIC="numeric"))
  #   results <- structure(list('estimates'=stime, 'loglik'=loglik,'niter'=nniter,'varmat'=varmat,'BIC'=round(BICCUB00,digits=8)), class = "cub")
  #   setMethod("LogLik",  function(cub) {
  #     cub@loglik
  #   })
  results<-list('estimates'=stime, 'loglik'=loglik,'niter'=nniter,'varmat'=varmat,'BIC'=round(BICCUB00,digits=8))
}
