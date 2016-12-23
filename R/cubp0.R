#' @title Main function for CUB models with covariates for the uncertainty component
#' @description Estimate and validate a CUB model for given ordinal responses, with covariates for explaining 
#' the feeling component via a logistic transform.
#' @aliases cubp0
#' @usage cubp0(m, ordinal, Y, maxiter, toler)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param Y Matrix of selected covariates for explaining the uncertainty component
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
#' @param toler Fixed error tolerance for final estimates 
#' @return An object of the class "CUB"
#' @import stats graphics
#' @references
#' Iannario M. and Piccolo D. (2010), A new statistical model for the analysis of customer satisfaction, 
#' \emph{Quality Technology and Quantity management}, \bold{7}(2) 149--168 \cr
#' Iannario M. and Piccolo D. (2012). CUB models: Statistical methods and empirical evidence, in:
#'  Kenett R. S. and Salini S. (eds.), \emph{Modern Analysis of Customer Surveys: with applications using R}, 
#'  J. Wiley and Sons, Chichester, 231--258
#' @keywords internal 

cubp0<-function(m,ordinal,Y,maxiter,toler){
  tt0<-proc.time()
  n<-length(ordinal)
  Y<-as.matrix(Y)
  
  if (ncol(Y)==1){
    Y<-as.numeric(Y)
  }
  
  p<-NCOL(Y)
  aver<-mean(ordinal); varcamp<-mean(ordinal^2)-aver^2;
  YY<-cbind(1,Y)
  ##################################################################
  serie<-1:m; freq<-tabulate(ordinal,nbins=m);
  inipaicsi<-inibest(m,freq)
  pai<-inipaicsi[1]; bet0<-log(pai/(1-pai));  
  betjj<- c(bet0,rep(0.1,p))                #betjj<-rep(0.1,p+1);
  csijj<-inipaicsi[2]
  ##############################################################
  loglikjj<-loglikcubp0(m,ordinal,Y,betjj,csijj)
  # ********************************************************************
  # ************* E-M algorithm for CUB(p,0) ***************************
  # ********************************************************************
  nniter<-1
  while(nniter<=maxiter){
    loglikold<-loglikjj
    bb<-probbit(m,csijj)
    vettn<-bb[ordinal]      # probbit for all ordinal (r_i,i=1,2,...,n)
    aai<- -1+ 1/(logis(Y,betjj)) #exp(-(YY%*%betjj));
    ttau<-1/(1+aai/(m*vettn))       # tau is a reserved word in R
    averpo<-sum(ordinal*ttau)/sum(ttau)
    ################################## maximize w.r.t. bet  ########
    bet<-betjj
    covar<-YY
    tauno<-ttau
    #nlmaxbet<-nlm(effe10,betjj,esterno10);   
    opmaxbet<-optim(bet,effe10,esterno10=cbind(tauno,covar))
    ################################################################         
    betjj<-opmaxbet$par
    # betjj<-nlmaxbet$estimate;        #updated bet estimates
    csijj<-(m-averpo)/(m-1)       #updated csi estimate
    #loglikjj<- -opmaxbet$value
    loglikjj<-loglikcubp0(m,ordinal,Y,betjj,csijj)
    
    #print(c(nniter,betjj,csijj,loglikjj)); #OPTIONAL PRINTING OF ITERATIONS
    testll<-abs(loglikjj-loglikold)
    if(testll<=toler) break else {loglikold<-loglikjj}
    nniter<-nniter+1
  }
  bet<-betjj;  csi<-csijj;  loglik<-loglikjj;
  ####################################################################
  AICCUBp0<- -2*loglik+2*(p+2)
  BICCUBp0<- -2*loglik+log(n)*(p+2)
  ####################################################################
  # Compute asymptotic standard errors of ML estimates
  ####################################################################
  varmat<-varcovcubp0(m,ordinal,Y,bet,csi) 
  nomi<-c(paste("beta",0:(length(bet)-1),sep="_"),"csi   ")
  stime<-c(bet,csi)
  nparam<-length(stime)
  #if(det(varmat)<=0) stop("Variance-covariance matrix NOT positive definite")
  if (isTRUE(varmat==matrix(NA,nrow=nparam,ncol=nparam))==TRUE){
    ddd<-cormat<-matrix(NA,nrow=nparam,ncol=nparam)
    ICOMP<-trvarmat<-NA
    errstd<-wald<-pval<-rep(NA,nparam) 
  } else {
    ddd<-diag(sqrt(1/diag(varmat)))
    cormat<-(ddd%*%varmat)%*%ddd
    trvarmat<-sum(diag(varmat))
    ICOMP<- -2*loglik + nparam*log(trvarmat/nparam) - log(det(varmat)) ## added
    errstd<-sqrt(diag(varmat));  wald<-stime/errstd;
    pval<-2*(1-pnorm(abs(wald)))
  }
  rownames(cormat)<-nomi;colnames(cormat)<-nomi; 
  
  durata<-proc.time()-tt0;durata<-durata[1];
  
  results<-list('estimates'=stime,'ordinal'=ordinal,'time'=durata,
                'loglik'=loglik,'niter'=nniter,'varmat'=varmat,
                'BIC'=BICCUBp0)
  #class(results)<-"cub"
  return(results)
}