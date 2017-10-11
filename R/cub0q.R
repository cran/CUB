#' @title Main function for CUB models with covariates for the feeling component
#' @description Function to estimate and validate a CUB model for given ordinal responses, with covariates for
#'  explaining the feeling component.
#' @aliases cub0q
#' @usage cub0q(m, ordinal, W, maxiter, toler)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param W Matrix of selected covariates for explaining the feeling component, not including intercept
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
#' @param toler Fixed error tolerance for final estimates 
#' @import stats graphics
#' @return An object of the class "CUB"
#' @references
#'  Piccolo D. and D'Elia A. (2008), A new approach for modelling consumers' preferences,
#'   \emph{Food Quality and Preference}, \bold{18}, 247--259 \cr
#' @references 
#' Iannario M. and Piccolo D. (2010), A new statistical model for the analysis of customer
#'  satisfaction, #' \emph{Quality Technology and Quantity management}, \bold{7}(2) 149--168 \cr
#' Iannario M. and Piccolo D. (2012), CUB models: Statistical methods and empirical evidence, in:
#' Kenett R. S. and Salini S. (eds.), \emph{Modern Analysis of Customer Surveys: with applications using R},
#'  J. Wiley and Sons, Chichester, 231--258. 
#' @keywords internal 
#' @examples 
#' #running donttest option since the proposed examples require a long time run for check 
#' \donttest{
#' data(relgoods)
#' m=10
#' ordinal=relgoods[,29]
#' gender=relgoods[,2]
#' data=na.omit(cbind(ordinal,gender))
#' ordinal=data[,1]
#' W=data[,2]
#' makeplot=TRUE
#' maxiter=500           
#' toler=1e-6 
#' cubfit=cub0q(m,ordinal,W, maxiter, toler, makeplot, summary=TRUE)
#' param=cubfit$estimates      # Final ML estimates
#' pai=param[1]                # Estimated uncertainty parameter
#' gama=param[2:length(param)] # Estimated coefficients for feeling covariates
#' maxlik=cubfit$loglik
#' varmat=cubfit$varmat
#' niter=cubfit$niter
#' BIC=cubfit$BIC
#' ###########################
#' data(univer)
#' m=7
#' global=univer[,12]
#' freqserv=univer[,2]
#' vercub0q=cub0q(m,global,W=freqserv,maxiter=300,toler=1e-4,makeplot=FALSE)
#' param=vercub0q$estimates      # Final ML estimates
#' pai=param[1]                  # Estimated uncertainty parameter
#' gama=param[2:length(param)]   # Estimated coefficients for feeling covariates
#' maxlik=vercub0q$loglik
#' varmat=vercub0q$varmat
#' niter=vercub0q$niter
#' BIC=vercub0q$BIC
#' }

cub0q<-function(m,ordinal,W,maxiter,toler){
  tt0<-proc.time()
  n<-length(ordinal)
  W<-as.matrix(W)
  if (ncol(W)==1){
    W<-as.numeric(W)
  }
  q<-NCOL(W)
  aver<-mean(ordinal)
  WW<-cbind(1,W)                   
  ##############################################################
  freq<-tabulate(ordinal,nbins=m); inipaicsi<-inibest(m,freq);  paijj<-inipaicsi[1]; 
  gamajj<-inibestgama(m,ordinal,W)
  ##############################################################
  loglikjj<-loglikcub0q(m,ordinal,W,paijj,gamajj)
  # ********************************************************************
  # ************* E-M algorithm for CUB(0,q) ***************************
  # ********************************************************************
  nniter<-1
  while(nniter<=maxiter){

    loglikold<-loglikjj
    vettn<-bitgama(m,ordinal,W,gamajj)
    ttau<-1/(1+(1-paijj)/(m*paijj*vettn)) 
    ################################# maximize w.r.t. gama  ########
    ordd<-ordinal;covar<-WW;
    gama<-gamajj
    optimgama<-optim(gama,effe01,esterno01=cbind(ttau,ordinal,WW),m=m) 
    ################################################################
    gamajj<-optimgama$par
    paijj<-sum(ttau)/n                    #updated pai estimate
 
    
    loglikjj<-loglikcub0q(m,ordinal,W,paijj,gamajj)## needed for nlm version
    # print(c(nniter,paijj,gamajj,loglikjj)); #OPTIONAL PRINTING OF ITERATIONS
    testll<-abs(loglikjj-loglikold)
    
   
    if(testll<=toler) break else {loglikold<-loglikjj}
    nniter<-nniter+1
  }
  pai<-paijj;  gama<-gamajj;  loglik<-loglikjj;
  ####################################################################
  AICCUB0q<- -2*loglik+2*(q+2)
  BICCUB0q<- -2*loglik+log(n)*(q+2)
  ####################################################################
  # Compute asymptotic standard errors of ML estimates
  ####################################################################
  varmat<-varcovcub0q(m,ordinal,W,pai,gama)
  nomi<-c("pai    ",paste("gamma",0:(length(gama)-1),sep="_"))
  stime<-c(pai,gama)
  nparam<-length(stime)
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
                'BIC'=BICCUB0q)
  # class(results)<-"cub"
  return(results)
}