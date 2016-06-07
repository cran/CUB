#' @title Main function for CUB models with covariates for both the uncertainty and the feeling components
#' @description Estimate and validate a CUB model for given ordinal responses, with covariates for explaining both the
#'  feeling and the uncertainty components by means of logistic transform.
#' @aliases cubpq
#' @usage cubpq(m, ordinal, Y, W, maxiter, toler, summary)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param Y Matrix of selected covariates for explaining the uncertainty component
#' @param W Matrix of selected covariates for explaining the feeling component
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
#' @param toler Fixed error tolerance for final estimates 
#' @param summary Logical: if TRUE, summary results of the fitting procedure are displayed on screen
#' @return An object of the class "CUB"
#' @import stats
#' @seealso \code{\link{varcovcubpq}}, \code{\link{loglikcubpq}}, \code{\link{inibestgama}}, \code{\link{CUB}}
#' @references
#' Piccolo D. and D'Elia A. (2008), A new approach for modelling consumers' preferences, \emph{Food Quality and Preference},
#' \bold{18}, 247--259 \cr
#' Iannario M. and Piccolo D. (2010), A new statistical model for the analysis of customer satisfaction, 
#' \emph{Quality Technology and Quantitative Management}, \bold{17}(2)  149--168
#' @keywords internal 

cubpq<-function(m,ordinal,Y,W,maxiter,toler,summary){
  tt0<-proc.time()
  n<-length(ordinal)
  p<-NCOL(Y)
  q<-NCOL(W)
  aver<-mean(ordinal)
  YY<-cbind(1,Y);     WW<-cbind(1,W);          
  #################################################################################
  freq<-tabulate(ordinal,nbins=m)
  inipaicsi<-inibest(m,freq); pai<-inipaicsi[1]; bet0<-log(pai/(1-pai)); betjj<-c(bet0,rep(0.1,p));
  gamajj<-inibestgama(m,ordinal,W)  ### Attention !!!
  #################################################################################
  loglikjj<-loglikcubpq(m,ordinal,Y,W,betjj,gamajj)
  # ********************************************************************
  # ************* E-M algorithm for CUB(p,q) ***************************
  # ********************************************************************
  nniter<-1
  while(nniter<=maxiter){
    loglikold<-loglikjj
    vettn<-bitgama(m,ordinal,W,gamajj)   
    aai<- -1+1/(logis(Y,betjj))   
    ttau<-1/(1+aai/(m*vettn))        
    ####################  maximize w.r.t. bet and gama    ############
    esterno10<-cbind(ttau,YY)
    esterno01<-cbind(ttau,ordinal,WW)
    bet<-betjj;  gama<-gamajj;
    betoptim<-optim(bet,effe10,esterno10=esterno10)
    gamaoptim<-optim(gama,effe01,esterno01=esterno01,m=m)
    ################################################################         
    betjj<-betoptim$par
    gamajj<-gamaoptim$par
    loglikjj<-loglikcubpq(m,ordinal,Y,W,betjj,gamajj)
    # print(c(nniter,betjj,gamajj,loglikjj)); #OPTIONAL PRINTING OF ITERATIONS
    testll<-abs(loglikjj-loglikold)
    if(testll<=toler) break else {loglikold<-loglikjj}
    nniter<-nniter+1
  }
  bet<-betjj;  gama<-gamajj;  loglik<-loglikjj;
  ####################################################################
  AICCUBpq<- -2*loglik+2*(p+q+2)
  BICCUBpq<- -2*loglik+log(n)*(p+q+2)
  ####################################################################
  # Compute asymptotic standard errors of ML estimates
  ####################################################################
  varmat<-varcovcubpq(m,ordinal,Y,W,bet,gama)
  #if(det(varmat)<=0) stop("Variance-covariance matrix NOT positive definite")
  nomi<-c(paste("beta",0:(length(bet)-1),sep="_"),paste("gamma",0:(length(gama)-1),sep="_"))
  stime<-round(c(bet,gama),5)
  nparam<-length(stime)
  if (isTRUE(varmat==matrix(NA,nrow=nparam,ncol=nparam))==TRUE){
    ddd<-cormat<-matrix(NA,nrow=nparam,ncol=nparam)
    ICOMP<-trvarmat<-NA
    errstd<-wald<-pval<-rep(NA,nparam)
  } else {
    ddd<-diag(sqrt(1/diag(varmat)))
    cormat<-(ddd%*%varmat)%*%ddd  
    trvarmat<-sum(diag(varmat))
    ICOMP<- -2*loglik + nparam*log(trvarmat/nparam) - log(det(varmat))
    errstd<-round(sqrt(diag(varmat)),5);  wald<-round(stime/errstd,5);
    pval<-round(2*(1-pnorm(abs(wald))),20)
  }
  rownames(cormat)<-nomi;colnames(cormat)<-nomi; 
  durata<-proc.time()-tt0;durata<-durata[1];
  ####################################################################
  # Print CUB(p,q) results of ML estimation  
  ####################################################################
  if (summary==TRUE){
    cat("\n")
    cat("==========================================================================","\n")
    cat("=======>>> C U B (p,q) model <<<=========   ML-estimates via E-M algorithm   ","\n")
    cat("==========================================================================","\n")
    cat(" Covariates for pai ==> p=",p,"     and     Covariates for csi ==> q=", q,"\n")
    cat("==========================================================================","\n")
    cat("*** m=", m,"  *** Sample size: n=", n,"   *** Iterations=",nniter,"Maxiter=",maxiter,"\n")
    cat("==========================================================================","\n")
    cat("parameters  ML-estimates  stand.errors    Wald-test      p-value ","\n")
    cat("==========================================================================","\n")
    for(i in 1:length(nomi)){
      cat(nomi[i],"     ",stime[i],"      ",errstd[i],"       ",wald[i],"      ",pval[i],"\n")
    }
    ####################################################################
    cat("==========================================================================","\n")
    cat("                         Parameters correlation matrix","\n") 
    print(round(cormat,5))
    ##############################################################################
    cat("==========================================================================","\n")
    cat("Log-lik(beta^,gamma^)=",round(loglik,digits=8),"\n")
    cat("Mean Log-likelihood  =",round(loglik/n,digits=8),"\n")
    cat("--------------------------------------------------------------------------","\n")
    cat("AIC-CUBpq           =",round(AICCUBpq,digits=8),"\n")
    cat("BIC-CUBpq           =",round(BICCUBpq,digits=8),"\n")
    cat("ICOMP-CUBpq         =",round(ICOMP,digits=8),"\n")
    cat("==========================================================================","\n")
    cat("Elapsed time              =",durata,"seconds","=====>>>",date(),"\n")
    cat("==========================================================================","\n")  
    
  }

  results<-list('estimates'=stime,'loglik'=loglik,'niter'=nniter,'varmat'=varmat, 'BIC'=round(BICCUBpq,digits=8))
}