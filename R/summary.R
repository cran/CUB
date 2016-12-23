#' @title S3 method: summary for class "GEM"
#' @description S3 method summary for objects of class \code{\link{GEM}}. 
#' @aliases summary.GEM
#' @method summary GEM
#' @param object An object of class \code{\link{GEM}}
#' @param correlation Logical: should the estimated correlation matrix be returned? Default is FALSE
#' @param ...  Other arguments
#' @export 
#' @return Extended summary results of the fitting procedure, including parameter estimates, their standard errors and
#' Wald statistics, maximized log-likelihood compared with that of the saturated model and of a Uniform sample.
#' AIC, BIC and ICOMP indeces are also displayed for model selection. Execution time and number of exectued iterations 
#' for the fitting procedure are aslo returned.
#' @import methods
#' @rdname summary.GEM
#' @keywords package
#' @examples 
#' data(relgoods)
#' attach(relgoods)
#' ordinal<-na.omit(MeetRelatives) 
#' model<-GEM(Formula(ordinal~0|0|0),family="cube") 
#' summary(model,correlation=TRUE,digits=4)
#'

# summary.GEM <- function(object, correlation=FALSE,  ...){
# 
#   if(!is.null(cl <- object$call)) {
#     cat("Call:\n")
#     dput(cl, control = NULL)
#   }
# 
#   ellipsis<-object$ellipsis
#   maxiter<-ellipsis$maxiter
#   family<-object$family
#   niter<-object$niter
# 
#   modello<-object$formula
#   #EFFE<-mod$Formula
#   data<-object$data
# 
#   m<-ellipsis$m
#   n<-length(object$ordinal)
# 
#   # if (missing(digits)){
#   #   digits<-5
#   # }
# 
#   stime<-round(object$estimates,5);
#   nomi<-parnames(object)
# 
# 
#   ordinal <- unclass(ordinal)
#   freq<-tabulate(ordinal,nbins=m)
#   varm<- vcov(object)    #as.matrix(object$varmat);
#   np<-length(stime)
#   if (isTRUE(varm==matrix(NA,nrow=np,ncol=np))==TRUE){
#     trvarmat<-ICOMP<-NA
#     errstd<-wald<-pval<-rep(NA,np)
#   } else {
#     trvarmat<-sum(diag(varm))
#     loglik<- as.numeric(logLik(object))
#     ICOMP<- -2*loglik + np*log(trvarmat/np) - log(det(varm))
#     errstd<-round(sqrt(diag(varm)),5);wald<-round(stime/errstd,5);
#     pval<-round(2*(1-pnorm(abs(wald))),20)
#   }
# 
#   AIC<- -2*loglik+2*(np)
#   BIC<- -2*loglik+log(n)*(np)
# 
#   llunif<- -n*log(m);
# 
#   nonzero<-which(freq!=0)
#   logsat <- -n*log(n)+sum((freq[nonzero])*log(freq[nonzero]))
#   devian<-2*(logsat-loglik)
# 
#   durata<-object$time
#   cat("\n")
#   cat("=======================================================================","\n")
#   cat("=====>>>", family," model    <<<=====   ML-estimates via E-M algorithm  ","\n")
#   cat("=======================================================================","\n")
#   cat(" m=", m," Sample size: n=", n," Iterations=",niter," Maxiter=",maxiter,"\n")
#   cat("=======================================================================","\n")
#  # cat("parameters  ML-estimates  stand.errors    estimates/stand.errors       ","\n")
#  #  cat("=======================================================================","\n")
#  #  estimates/stand.errors
#   StdErr<-errstd
#   Wald<-wald
#  # matout<-cbind(nomi,stime,errstd,Wald)
#  #   dimnames(matout)<-list(rep("",length(nomi)),c("Parameters",  "ML-estimates" , "Std. err.", "Est./Std.err (Wald test)"))
# 
#   matout<-cbind(coef(object),StdErr,Wald)
#   print(matout)
#     # for(i in 1:np){
#   #   cat(nomi[i],"      ",stime[i],"        ",errstd[i],"  ",wald[i],"      ","\n")
#   # }
#   cat("=======================================================================","\n")
#   if (family=="CUB" & !is.null(ellipsis$shelter)){
#     covshe<-model.matrix(modello,data=data,rhs=3)
#     if (ncol(covshe)==0){
#       X<-NULL
#     } else {
#       X<-covshe[,-1]
#     }
#     if (is.null(X)){
#       pai1<-stime[1];pai2<-stime[2];csi<-stime[3]
#       delta<-1-pai1-pai2
#       paistar<-pai1/(pai1+pai2)
#       stime2<-c(paistar,csi,delta)
#       nomi2<-c("paistar","csi","delta")
#       vv<-varm
#       esdelta<-sqrt(vv[1,1]+vv[2,2]+2*vv[1,2])
#       espaistar<-paistar*(1-paistar)*sqrt(vv[1,1]/(pai1^2) -2*vv[1,2]/(pai2*pai1) + vv[2,2]/pai2^2   )
#       errstd2<-c(espaistar,errstd[3],esdelta)
#       wald2<-round(stime2/errstd2,5);
#       ErrStd<-as.numeric(errstd2)
# 
#       mat2<-as.matrix(stime2)
#       dimnames(mat2)<-list(nomi2,"")
#       Wald<-wald2
#       matout2<-cbind(mat2,ErrStd,wald2)
#       #matout2<-cbind(nomi2,stime2,errstd2,wald2)
#     #  dimnames(matout2)<-list(rep("",length(nomi2)),c("Parameters",  "ML-estimates" , "Std. err.", "Est./Std.err (Wald test)"))
#       cat("Alternative parameterization","\n")
#       print(matout2)
#       # for(i in 1:np){
#       #   cat(nomi2[i],"      ",stime2[i],"         ",errstd2[i],"       ",wald2[i],"      ","\n")
#       # }
#       cat("=======================================================================","\n")
#     }
#   }
# 
#   if (correlation==TRUE){
#     cat("Parameters Correlation matrix","\n")
#     print(cormat(object))
#     cat("=======================================================================","\n")
#   }
# 
#   cat("Log-lik            =",round(loglik,digits=7),"\n")
#   cat("Mean Log-likelihood=",round(loglik/n,digits=7),"\n")
#   cat("Log-lik(UNIFORM)   =",round(llunif,digits=7),"\n")
#   cat("Log-lik(saturated) =",round(logsat,digits=7),"\n")
#   cat("Deviance           =",round(devian,digits=7),"\n")
#   cat("-----------------------------------------------------------------------","\n")
#   #cat("Log-lik(Shifted-BINOMIAL) =",round(llsb,digits=8),"\n")
#   #cat("-----------------------------------------------------------------------","\n")
#   cat("AIC       =",round(AIC,digits=7),"\n")
#   cat("BIC       =",round(BIC,digits=7),"\n")
#   cat("ICOMP     =",round(ICOMP,digits=7),"\n")
#   cat("=======================================================================","\n")
#   cat("Elapsed time=",object$time,"seconds","=====>>>",date(),"\n")
#   cat("=======================================================================","\n")
# 
# 
# }




################################################################
###############################################################

#summary <- function(object,...) UseMethod("summary", object)

# digits=options()$digits 

summary.GEM <- function(object, correlation=FALSE, ...){

  
  arguments<-list(...)
  
  digits<-arguments$digits
  
  if (is.null(digits)){
    digits<-options()$digits
  }
  
  ellipsis<-object$ellipsis
  m<-ellipsis$m
  n<-length(object$ordinal)

  output<-list()

  stime<-object$estimates
  ordinal <- unclass(object$ordinal)
  freq<-tabulate(ordinal,nbins=m)
  varm<- vcov(object)    #as.matrix(object$varmat);
  np<-length(stime)
  if (isTRUE(varm==matrix(NA,nrow=np,ncol=np))==TRUE){
    trvarmat<-output$ICOMP<-NA
    output$errstd<-output$wald<-output$pval<-rep(NA,np)
  } else {
    trvarmat<-sum(diag(varm))
    output$loglik<- as.numeric(logLik(object))
    output$ICOMP<- -2*output$loglik + np*log(trvarmat/np) - log(det(varm))
    output$errstd<-sqrt(diag(varm));
    output$wald<-stime/output$errstd;
    output$pval<-round(2*(1-pnorm(abs(output$wald))),digits)
  }

  output$loglik<-logLik(object)
  output$AIC<- -2*logLik(object)+2*(np)
  output$BIC<- -2*logLik(object)+log(n)*(np)

  output$llunif<- -n*log(m);

  nonzero<-which(freq!=0)
  output$logsat <- -n*log(n)+sum((freq[nonzero])*log(freq[nonzero]))
  output$devian<-2*(output$logsat-object$loglik)
  output$object<-object
  output$n<-n

  output$cormat<-NULL
  if (correlation==TRUE){
    output$cormat<- cormat(object)
  }

  class(output)<-"summary.GEM"
  
  
  
  print.summary.GEM <- function(x,...){
    
    if(!is.null(cl <- x$call)) {
      cat("Call:\n")
      dput(cl, control = NULL)
    }
    
    ellipsis<-x$ellipsis
    object<-x$object
    maxiter<-object$ellipsis$maxiter
    family<-object$family
    niter<-object$niter
    m<-object$ellipsis$m
    
    modello<-object$formula
    data<-object$data
    n<-x$n
 #   cat("\n")
    cat("=======================================================================","\n")
    cat("=====>>>", family," model    <<<=====   ML-estimates via E-M algorithm  ","\n")
    cat("=======================================================================","\n")
    cat(" m=", m," Sample size: n=",n," Iterations=", niter," Maxiter=",maxiter,"\n")
    cat("=======================================================================","\n")
    # cat("parameters  ML-estimates  stand.errors    estimates/stand.errors       ","\n")
    #  cat("=======================================================================","\n")
    #  estimates/stand.errors
    StdErr<-x$errstd
    Wald<-x$wald
    # matout<-cbind(nomi,stime,errstd,Wald)
    #   dimnames(matout)<-list(rep("",length(nomi)),c("Parameters",  "ML-estimates" , "Std. err.", "Est./Std.err (Wald test)"))
    
    matout<-cbind(coef(object),StdErr,Wald)
    colnames(matout)<-c("Estimates","StdErr","Wald")
    print(matout,digits=digits)
    # for(i in 1:np){
    #   cat(nomi[i],"      ",stime[i],"        ",errstd[i],"  ",wald[i],"      ","\n")
    # }
    cat("=======================================================================","\n")
    if (family=="CUB" & !is.null(object$ellipsis$shelter)){
      covshe<-model.matrix(modello,data=data,rhs=3)
      if (ncol(covshe)==0){
        X<-NULL
      } else {
        X<-covshe[,-1]
      }
      if (is.null(X)){
        #stime<-as.numeric(coef(object))
        pai1<-stime[1];pai2<-stime[2];csi<-stime[3]
        delta<-1-pai1-pai2
        paistar<-pai1/(pai1+pai2)
        stime2<-c(paistar,csi,delta)
        nomi2<-c("paistar","csi","delta")
        vv<-vcov(object)
        esdelta<-sqrt(vv[1,1]+vv[2,2]+2*vv[1,2])
        espaistar<-paistar*(1-paistar)*sqrt(vv[1,1]/(pai1^2) -2*vv[1,2]/(pai2*pai1) + vv[2,2]/pai2^2)
        errstd2<-c(espaistar,StdErr[3],esdelta)
        wald2<-stime2/errstd2
        ErrStd<-as.numeric(errstd2)
        
        mat2<-as.matrix(stime2)
        dimnames(mat2)<-list(nomi2,"")
        Wald<-wald2
        matout2<-cbind(mat2,ErrStd,wald2)
        colnames(matout2)<-c("Estimates","StdErr","Wald")
        #matout2<-cbind(nomi2,stime2,errstd2,wald2)
        #  dimnames(matout2)<-list(rep("",length(nomi2)),c("Parameters",  "ML-estimates" , "Std. err.", "Est./Std.err (Wald test)"))
        cat("Alternative parameterization","\n")
        print(matout2,digits=digits)
        # for(i in 1:np){
        #   cat(nomi2[i],"      ",stime2[i],"         ",errstd2[i],"       ",wald2[i],"      ","\n")
        # }
        cat("=======================================================================","\n")
      }
    }
    
    if (!is.null(x$cormat)){
      cat("Parameters Correlation matrix","\n")
      print(x$cormat)
      cat("=======================================================================","\n")
    }
    loglik<-logLik(object)
    cat("Log-lik            =",round(loglik,digits=digits),"\n")
    cat("Mean Log-likelihood=",round(loglik/n,digits=digits),"\n")
    cat("Log-lik(UNIFORM)   =",round(x$llunif,digits=digits),"\n")
    cat("Log-lik(saturated) =",round(x$logsat,digits=digits),"\n")
    cat("Deviance           =",round(x$devian,digits=digits),"\n")
    cat("-----------------------------------------------------------------------","\n")
    #cat("Log-lik(Shifted-BINOMIAL) =",round(llsb,digits=8),"\n")
    #cat("-----------------------------------------------------------------------","\n")
    cat("AIC       =",round(x$AIC,digits=digits),"\n")
    cat("BIC       =",round(x$BIC,digits=digits),"\n")
    cat("ICOMP     =",round(x$ICOMP,digits=digits),"\n")
    cat("=======================================================================","\n")
    cat("Elapsed time=",object$time,"seconds","=====>>>",date(),"\n")
    cat("=======================================================================","\n")
    
    
  }
  
  
  
  print(output)
  #output
   #invisible(output)

}

