#' @title Main function for CUBE models with covariates
#' @description Function to estimate and validate a CUBE model with 
#' explicative covariates for all the three parameters.
#' @aliases cubecov
#' @keywords internal
#' @usage cubecov(m, ordinal, Y, W, Z, starting, maxiter, toler)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param Y Matrix of selected covariates for explaining the uncertainty component
#' @param W Matrix of selected covariates for explaining the feeling component
#' @param Z Matrix of selected covariates for explaining the overdispersion component
#' @param starting Vector of initial parameters estimates to start the optimization algorithm 
#' (it has length NCOL(Y) + NCOL(W) + NCOL(Z) + 3 to account for intercept terms 
#' for all the three components
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
#' @param toler Fixed error tolerance for final estimates 
#' @return An object of the class "CUBE"
#' @import stats
#' @references
#' Piccolo, D. (2014). Inferential issues on CUBE models with covariates,
#'  \emph{Communications in Statistics - Theory and Methods}, \bold{44}, 
#'  DOI: 10.1080/03610926.2013.821487

cubecov<-function(m,ordinal,Y,W,Z,starting,maxiter,toler){
  tt0<-proc.time()
  n<-length(ordinal)
  Y<-as.matrix(Y); W<-as.matrix(W); Z<-as.matrix(Z);
  p<-NCOL(Y)
  q<-NCOL(W)
  v<-NCOL(Z)
  aver<-mean(ordinal)
  if (ncol(W)==1){
    W<-as.numeric(W)
  }
  if (ncol(Y)==1){
    Y<-as.numeric(Y)
  }
  if (ncol(Z)==1){
    Z<-as.numeric(Z)
  }
  
  
  YY<-cbind(1,Y)   
  WW<-cbind(1,W)       
  ZZ<-cbind(1,Z)    
  # **********************************
  # *** E-M algorithm for CUBECOV  ***
  # **********************************
  #################################################################
  # 00.# Initial values of parameters ............. Attention !!! #
  #################################################################
  betjj<-starting[1:(p+1)]
  gamajj<-starting[(p+2):(p+q+2)]
  alphajj<-starting[(p+q+3):(p+q+v+3)]
  ################################
  # * * * Iterative Loop    * * *#
  ################################
  nniter<-1
  while(nniter<=maxiter){
    #################################################################
    # 01.# Initial values of vectors i=1,2,..,n
    #################################################################
    paijj<-logis(Y,betjj); csijj<-logis(W,gamajj) 
    phijj<-1/(1/logis(Z,alphajj)-1)  #****
    #################################################################
    # 02.# Computation of beta-binomial distribution i=1,2,..,n
    #################################################################
    betabin<-betabinomial(m,factor(ordinal,ordered=TRUE),csijj,phijj)
    #################################################################
    # 03.# Computation of CUBE probability distribution i=1,2,..,n
    #################################################################
    probi<-paijj*(betabin-1/m)+1/m
    likold<-sum(log(probi))
    #################################################################
    # 4.# Computation of conditional probability i=1,2,..,n
    #################################################################
    taui<-1-(1-paijj)/(m*probi)
    #################################################################
    # 5. Unify parameter vectors
    param<-c(gamajj,alphajj)
    #################################################################
    # 6.# Maximization of Q_1(beta) and Q_2(param)
    #################################################################
    ### maximize w.r.t. bet and gama #########
    esterno1<-cbind(taui,YY) 
    covar<-esterno1[,2:NCOL(esterno1)]
    esterno2<-cbind(taui,ordinal,W,Z) 
    bet<-betjj
    optimbet<-optim(bet,Quno,esterno1=esterno1,gr=NULL) #added gr  
    optimparam<-optim(param,Qdue,esterno2=esterno2,q=q,m=m,gr=NULL)
    #################################################################
    # 7.# Computation of updated estimates and log-likelihood
    #################################################################
    betjj<-optimbet$par                #updated bet estimates
    paramjj<-optimparam$par
    gamajj<-paramjj[1:(q+1)]           #updated gama estimates
    alphajj<-paramjj[(q+2):(q+v+2)]    #updated alpha estimates
    ### updated log-likelihood
    liknew<-loglikcubecov(m,ordinal,Y,W,Z,betjj,gamajj,alphajj)
    #################################################################
    # 8.# Checking improvement of updated log-likelihood
    #################################################################
    # print(c(nniter,betjj,gamajj,loglikjj)); #OPTIONAL PRINTING OF ITERATIONS
    testloglik<-abs(liknew-likold)
    #print(nniter);##added
    #print(round(c(liknew,likold,testloglik),digits=7))#added
    if(testloglik<=toler) break else {likold<-liknew}
    nniter<-nniter+1
  }
  #################################################################
  # 8.# Final ML estimates and maximized log-likelihood
  #################################################################
  bet<-betjj;  gama<-gamajj;  alpha<-alphajj;
  loglik<-liknew
  ###### End of E-M algorithm for CUBE ***********************************************
  paramest<-c(bet,gama,alpha)
  nparam<- length(paramest) ###p+q+v+3;
  AICCUBE<- -2*loglik+2*nparam
  BICCUBE<- -2*loglik+log(n)*nparam
  
  ############################################################
  # Compute asymptotic standard errors of ML CUBE estimates ##
  ############################################################
  
  varmat<-varcovcubecov(m,ordinal,Y,W,Z,bet,gama,alpha)
  #if(det(varmat)<=0) stop("Variance-Covariance matrix NOT positive definite")
  # nomi<-c(paste("beta",0:(length(bet)-1),sep="_"),
  #         paste("gamma",0:(length(gama)-1),sep="_"),
  #         paste("alpha",0:(length(alpha)-1),sep="_"))
  stime<-paramest
  
  if (isTRUE(varmat==matrix(NA,nrow=nparam,ncol=nparam))==TRUE){
    ddd<-cormat<-matrix(NA,nrow=nparam,ncol=nparam)
    trvarmat<-ICOMP<-NA
    errstd<-wald<-pval<-rep(NA,nparam)
  } else {
    ddd<-diag(sqrt(1/diag(varmat)))
    cormat<-(ddd%*%varmat)%*%ddd
    trvarmat<-sum(diag(varmat))
    ICOMP<- -2*loglik + nparam*log(trvarmat/nparam) - log(det(varmat))
    errstd<-sqrt(diag(varmat));  wald<-stime/errstd;
    pval<-2*(1-pnorm(abs(wald)))
  }
  # rownames(cormat)<-nomi;colnames(cormat)<-nomi; 
  ##################################################
  # Print CUBE-covariates results of ML estimation #
  durata<-proc.time()-tt0;durata<-durata[1];
  ##################################################
  # if (summary==TRUE){
  
  #cat("=======================================================================","\n")
  #cat("Convergence code =",optimparam$convergence,"\n")
  results<-list('estimates'=stime, 'loglik'=loglik, 'niter'= nniter,
                'varmat'=varmat,'BIC'=BICCUBE,'time'=durata)
  
}