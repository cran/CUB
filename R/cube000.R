#' @title Main function for CUBE models without covariates
#' @description Estimate and validate a CUBE model without covariates.
#' @aliases cube000
#' @usage cube000(m, ordinal, starting, maxiter, toler, expinform)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param starting Vector of initial estimates to start the optimization algorithm, 
#' whose length equals the number of parameters of the model
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
#' @param toler Fixed error tolerance for final estimates 
#' @param expinform Logical: if TRUE, the function returns the expected variance-covariance matrix
#' @return An object of the class "CUBE"
#' @import stats
#' @references
#'  Iannario, M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data,
#'   \emph{Communications in Statistics - Theory and Methods}, \bold{43}, 771--786 \cr
#' Iannario, M. (2015). Detecting latent components in ordinal data with overdispersion by means 
#' of a mixture distribution, \emph{Quality & Quantity}, \bold{49}, 977--987
#' @keywords internal #models

cube000<-function(m,ordinal,starting,maxiter,
                  toler,expinform){ #default for expinform = FALSE
  tt0<-proc.time()
  freq<-tabulate(ordinal,nbins=m); n<-sum(freq); 
  aver<-mean(ordinal); varcamp<-mean(ordinal^2)-aver^2;
  ########################################################
  #(00)# initial estimates, not efficient:
  #starting<-inibestcube(m,ordinal)
  pai<-starting[1]; csi<-starting[2]; phi<-starting[3];
  #(0)# log-lik
  loglik<-loglikcube(m,freq,pai,csi,phi)
  
  # ********************************************************************
  # *************   E-M algorithm for CUBE     *************************
  # ********************************************************************
  nniter<-1
  while(nniter<=maxiter){
    likold<-loglik
    
    #(1)# betar
    bb<-betar(m,csi,phi)
    aa<-(1-pai)/(m*pai*bb)
    
    #(2)# taunor
    tauno<-1/(1+aa)
    
    #(3)# pai(k+1)
    pai<-sum(freq*tauno)/n # updated pai estimate
    
    paravecjj<-c(csi,phi)
    
    #(4)# Q(k+1)
    dati<-cbind(tauno,freq)
    ################ EFFECUBE is Q(csi,phi) ###########################
    #(5)# (csi(k+1),phi(k+1))
    ################################## maximize w.r.t. paravec  ########
    paravec<-paravecjj
    optimestim<-optim(paravec,effecube,dati=dati,m=m,method = "L-BFGS-B",lower=c(0.01,0.01),upper=c(0.99,0.3))    # print(nlmaxg)
    ################################################################         
    
    #(6)# theta(k+1)
    paravecjj<-optimestim$par   # updated paravec estimates
    csi<-paravecjj[1];   phi<-paravecjj[2];
    ##########################################
    if(pai<0.001){pai<-0.001; nniter<-maxiter-1}
    #if(csi<0.001){csi<-0.001; nniter<-maxiter-1}
    #if(phi<0.001){phi<-0.001; nniter<-maxiter-1}
    if(pai>0.999){pai<-0.99}         ### to avoid division by 0 !!!
    #if(csi>0.999){csi<-0.99}         ### to avoid division by 0 !!!
    ###################################### print(c(nniter,pai,csi,phi));
    
    #(7)# elle(theta(k+1))
    liknew<-loglikcube(m,freq,pai,csi,phi)
    
    #(8)# test
    testll<-abs(liknew-likold)           # OPTIONAL printing: print(testll); 
    # OPTIONAL printing: print(cbind(nniter,testll,pai,csi,phi));
    if(testll<=toler) break else {loglik<-liknew} # OPTIONAL printing: print(loglik);
    nniter<-nniter+1
  }
  loglik<-liknew
  ###### End of E-M algorithm for CUBE ***********************************************
  AICCUBE<- -2*loglik+2*(3)
  BICCUBE<- -2*loglik+log(n)*(3)
  # ********************************************************
  # Compute ML var-cov matrix and print result for CUBE
  # ********************************************************
  
  if(expinform==TRUE){
    varmat<-varcovcubeexp(m,pai,csi,phi,n)
  }
  else{
    varmat<-varcovcubeobs(m,pai,csi,phi,freq)
  }
  
  # nomi<-rbind("pai  ","csi  ","phi  ")
  stime<-c(pai,csi,phi)
  nparam<-length(stime)
  if (isTRUE(varmat==matrix(NA,nrow=nparam,ncol=nparam))==TRUE){
    ddd<-matrix(NA,nrow=nparam,ncol=nparam)
    trvarmat<-ICOMP<-NA
    errstd<-wald<-pval<-rep(NA,nparam)
  } else {
    trvarmat<-sum(diag(varmat))
    ICOMP<- -2*loglik + nparam*log(trvarmat/nparam) - log(det(varmat)) 
    errstd<-sqrt(diag(varmat))
    wald<-stime/errstd
    pval<-2*(1-pnorm(abs(wald)))
    ddd<-diag(sqrt(1/diag(varmat)))
  }
  ####################################################################
  # Print CUBE results of ML estimation  
  ####################################################################
  ### Log-likelihood comparisons ##############################################
  llunif<- -n*log(m) 
  csisb<-(m-aver)/(m-1)
  llsb<-loglikcub00(m,freq,1,csisb)
  nonzero<-which(freq!=0)
  logsat<- -n*log(n)+sum((freq[nonzero])*log(freq[nonzero]))
  devian<-2*(logsat-loglik)
  theorpr<-probcube(m,pai,csi,phi)
  dissimcube<-dissim(theorpr,freq/n)
  pearson<-((freq-n*theorpr))/sqrt(n*theorpr)
  X2<-sum(pearson^2)
  relares<-(freq/n-theorpr)/theorpr
  FF2<-1-dissimcube
  LL2<-1/(1+mean((freq/(n*theorpr)-1)^2))
  II2<-(loglik-llunif)/(logsat-llunif)
  stampa<-cbind(1:m,freq/n,theorpr,pearson,relares)
  durata<-proc.time()-tt0;durata<-durata[1];
  
  results<-list('estimates'= stime, 'loglik'= loglik,
                'niter'= nniter, 'varmat'= varmat,'BIC'=BICCUBE,'time'=durata)
  
}