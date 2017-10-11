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
#' model<-GEM(Formula(MeetRelatives~0|0|0),family="cube",data=relgoods) 
#' summary(model,correlation=TRUE,digits=4)
#'


################################################################
###############################################################

#summary <- function(object,...) UseMethod("summary", object)

# digits=options()$digits 

summary.GEM <- function(object, correlation=FALSE, ...){

  flagcov<-0
  arguments<-list(...)
  
  digits<-arguments$digits
  
  if (is.null(digits)){
    digits<-options()$digits
  }
  
  ellipsis<-object$ellipsis
  m<-ellipsis[['m']]
  n<-length(object$ordinal)

  output<-list()

  stime<-object$estimates
 ordinal<-object$ordinal
  
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
  output$ellipsis<-ellipsis
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
 StdErr<-output$errstd
 Wald<-output$wald
 matout<-cbind(stime,StdErr,Wald)
 
 colnames(matout)<-c("Estimates","StdErr","Wald")
 rownames(matout)<-parnames(object)
 output$results<-matout
 
  class(output)<-"summary.GEM"
  
  
  
 print.summary.GEM <- function(x,...){
   
   if(!is.null(cl <- x$call)) {
     cat("Call:\n")
     dput(cl, control = NULL)
   }
   
   ellipsis<-x$ellipsis
   object<-x$object
   maxiter<-object$ellipsis[['maxiter']]
   family<-object$family
   niter<-object$niter
   m<-object$ellipsis[['m']]
   stime<-object$estimates
   
   
   modello<-object$formula
   data<-ellipsis$data
   
   mf<-model.frame(modello,data=data,na.action=na.omit)
   
   n<-x$n
   cat("=======================================================================","\n")
   cat("=====>>>", family," model    <<<=====   ML-estimates via E-M algorithm  ","\n")
   cat("=======================================================================","\n")
   cat(" m=", m," Sample size: n=",n," Iterations=", niter," Maxiter=",maxiter,"\n")
   cat("=======================================================================","\n")
   
   StdErr<-x$errstd
   Wald<-x$wald
   
   
   data<-object$data
   
   listanomi<-parnames(object)
   
   if ( family == "CUB"){
     covpai<-model.matrix(modello,data=mf,rhs=1)
     covcsi<-model.matrix(modello,data=mf,rhs=2)
     covshe<-model.matrix(modello,data=mf,rhs=3)
     
     if (ncol(covpai)!=0 | ncol(covcsi)!=0 | ncol(covshe)!=0){
       flagcov<-1
     }
     
     if (ncol(covpai)==0){
       Y<-NULL
     } else {
       if (NCOL(covpai)==2){
         Y<-as.matrix(covpai[,-1])
         colnames(Y)<-colnames(covpai)[2]
       } else {
         Y<-covpai[,-1]
       }
     }
     if (ncol(covcsi)==0){
       W<-NULL
     } else {
       if (NCOL(covcsi)==2){
         W<-as.matrix(covcsi[,-1])
         colnames(W)<-colnames(covcsi)[2]
       } else {
         W<-covcsi[,-1]
       }
     }
     
     if (ncol(covshe)==0){
       X<-NULL
     } else {
       X<-covshe[,-1]
     }
     
     if (!is.null(X) & !is.null(Y) & !is.null(W) & !is.null(object$ellipsis$shelter)){
       Y<-as.matrix(Y); W<-as.matrix(W); X<-as.matrix(X);
       p<-NCOL(Y);
       q<-NCOL(W); 
       s<-NCOL(X); 
       
       mat1<-cbind(stime[1:(p+1)],StdErr[1:(p+1)],Wald[1:(p+1)])
       colnames(mat1)<-c("Estimates","StdErr","Wald")
       rownames(mat1)<-listanomi[1:(p+1)]
       x$uncertainty<-mat1
       
       mat2<-cbind(stime[(p+2):(p+q+2)],StdErr[(p+2):(p+q+2)],Wald[(p+2):(p+q+2)])
       x$feeling<-mat2
       
       mat3<-cbind(stime[(p+q+3):(p+q+s+3)],StdErr[(p+q+3):(p+q+s+3)],Wald[(p+q+3):(p+q+s+3)])
       x$shelter<-mat3
       
       cat("Uncertainty                                           ", "\n")
       
       print(mat1,digits=digits)
       cat("=======================================================================","\n")
       cat("Feeling                                           ", "\n")
       colnames(mat2)<-c("Estimates","StdErr","Wald")
       rownames(mat2)<-listanomi[(p+2):(p+q+2)]
       
       print(mat2,digits=digits)
       cat("=======================================================================","\n")
       
       cat("Shelter effect                                ", "\n")
       colnames(mat3)<-c("Estimates","StdErr","Wald")
       rownames(mat3)<-listanomi[(p+q+3):(p+q+s+3)]
       print(mat3,digits=digits)
       
     } else if (is.null(object$ellipsis$shelter) & is.null(X) & !is.null(Y) & !is.null(W)){
       Y<-as.matrix(Y); W<-as.matrix(W);
       p<-NCOL(Y);
       q<-NCOL(W); 
       
       mat1<-cbind(stime[1:(p+1)],StdErr[1:(p+1)],Wald[1:(p+1)])
       colnames(mat1)<-c("Estimates","StdErr","Wald")
       rownames(mat1)<-listanomi[1:(p+1)]
       x$uncertainty<-mat1
       
       mat2<-cbind(stime[(p+2):(p+q+2)],StdErr[(p+2):(p+q+2)],Wald[(p+2):(p+q+2)])
       cat("Uncertainty                                           ", "\n")
       
       print(mat1,digits=digits)
       colnames(mat2)<-c("Estimates","StdErr","Wald")
       rownames(mat2)<-listanomi[(p+2):(p+q+2)]
       
       x$uncertainty<-mat1
       x$feeling<-mat2
       
       cat("=======================================================================","\n")
       
       cat("Feeling                                           ", "\n")
       print(mat2,digits=digits)
       
     } else if (is.null(object$ellipsis$shelter) & is.null(X) & is.null(Y) & !is.null(W)){
       W<-as.matrix(W);
       
       q<-NCOL(W); 
       
       mat1<-cbind(stime[1],StdErr[1],Wald[1])
       colnames(mat1)<-c("Estimates","StdErr","Wald")
       rownames(mat1)<-listanomi[1]
         
       mat2<-cbind(stime[2:(q+2)],StdErr[2:(q+2)],Wald[2:(q+2)])
       colnames(mat2)<-c("Estimates","StdErr","Wald")
       rownames(mat2)<-listanomi[2:(q+2)]
       cat("Uncertainty                                           ", "\n")
       
       print(mat1,digits=digits)
       cat("=======================================================================","\n")
       
       cat("Feeling                                           ", "\n")
       print(mat2,digits=digits)
       x$uncertainty<-mat1
       x$feeling<-mat2
     } else if (is.null(object$ellipsis$shelter) & is.null(X) & !is.null(Y) & is.null(W)){
       Y<-as.matrix(Y);
       
       p<-NCOL(Y); 
       
       mat1<-cbind(stime[1:(p+1)],StdErr[1:(p+1)],Wald[1:(p+1)])
       colnames(mat1)<-c("Estimates","StdErr","Wald")
       rownames(mat1)<-listanomi[1:(p+1)]
       
       mat2<-cbind(stime[(p+2)],StdErr[p+2],Wald[p+2])
       colnames(mat2)<-c("Estimates","StdErr","Wald")
       rownames(mat2)<-listanomi[p+2]
       cat("Uncertainty                                           ", "\n")
       
       print(mat1,digits=digits)
       cat("=======================================================================","\n")
       
       cat("Feeling                                           ", "\n")
       print(mat2,digits=digits)
       x$uncertainty<-mat1
       x$feeling<-mat2
     } else if (is.null(object$ellipsis$shelter) & is.null(X) & is.null(Y) & is.null(W)) {
       
       mat1<-cbind(stime[1],StdErr[1],Wald[1])
       colnames(mat1)<-c("Estimates","StdErr","Wald")
       
       mat2<-cbind(stime[2],StdErr[2],Wald[2])
       colnames(mat2)<-c("Estimates","StdErr","Wald")
       rownames(mat2)<-listanomi[2]
       x$uncertainty<-mat1
       x$feeling<-mat2
       
       cat("Uncertainty                                           ", "\n")
       
       print(mat1,digits=digits)
       cat("=======================================================================","\n")
       
       cat("Feeling                                           ", "\n")
       print(mat2,digits=digits)
       
       
     }
   }
   
   
   if (family == "CUBE"){
    
     covpai<-model.matrix(modello,data=mf,rhs=1)
     covcsi<-model.matrix(modello,data=mf,rhs=2)
     covphi<-model.matrix(modello,data=mf,rhs=3)
     
     if (ncol(covpai)!=0 | ncol(covcsi)!=0 | ncol(covphi)!=0){
       flagcov<-1
     }
     
     if (ncol(covpai)==0){
       Y<-NULL
     } else {
       Y<-covpai[,-1]
     }
     if (ncol(covcsi)==0){
       W<-NULL
     } else {
       W<-covcsi[,-1]
     }
     if (ncol(covphi)==0){
       Z<-NULL
     } else {
       Z<-covphi[,-1]
     }
     
     if (is.null(Y)& is.null(W) & is.null(Z)){
       mat1<-cbind(stime[1],StdErr[1],Wald[1])
       colnames(mat1)<-c("Estimates","StdErr","Wald")
       rownames(mat1)<-listanomi[1]
       mat2<-cbind(stime[2],StdErr[2],Wald[2])
       colnames(mat2)<-c("Estimates","StdErr","Wald")
       rownames(mat2)<-listanomi[2]
       
       cat("Uncertainty                                           ", "\n")
       
       print(mat1,digits=digits)
       cat("=======================================================================","\n")
       
       cat("Feeling                                           ", "\n")
       print(mat2,digits=digits)
       cat("=======================================================================","\n")
       
       cat("Overdispersion                                           ", "\n")
       
       mat3<-cbind(stime[3],StdErr[3],Wald[3])
       colnames(mat3)<-c("Estimates","StdErr","Wald")
       rownames(mat3)<-listanomi[3]
       print(mat3,digits=digits)
       
       x$uncertainty<-mat1
       x$feeling<-mat2
       x$overdispersion<-mat3
       
     } else if (is.null(Y)& !is.null(W) & is.null(Z)){
       
       q<-NCOL(W)
       
       mat1<-cbind(stime[1],StdErr[1],Wald[1])
       colnames(mat1)<-c("Estimates","StdErr","Wald")
       rownames(mat1)<-listanomi[1]
       
       mat2<-cbind(stime[2:(q+2)],StdErr[2:(q+2)],Wald[2:(q+2)])
       colnames(mat2)<-c("Estimates","StdErr","Wald")
       rownames(mat2)<-listanomi[2:(q+2)]
       
       cat("Uncertainty                                           ", "\n")
       print(mat1,digits=digits)
       
       cat("=======================================================================","\n")
       
       cat("Feeling                                           ", "\n")
       print(mat2,digits=digits)
       cat("=======================================================================","\n")
       
       cat("Overdispersion                                           ", "\n")
       
       mat3<-cbind(stime[q+3],StdErr[q+3],Wald[q+3])
       colnames(mat3)<-c("Estimates","StdErr","Wald")
       rownames(mat3)<-listanomi[q+3]
       
       print(mat3,digits=digits)
       x$uncertainty<-mat1
       x$feeling<-mat2
       x$overdispersion<-mat3
       
     }
     else if (!is.null(Y)& !is.null(W) & !is.null(Z)){
       p<-NCOL(Y)
       q<-NCOL(W)
       s<-NCOL(Z)
       
       mat1<-cbind(stime[1:(p+1)],StdErr[1:(p+1)],Wald[1:(p+1)])
       mat2<-cbind(stime[(2+p):(q+p+2)],StdErr[(2+p):(q+p+2)],Wald[(2+p):(q+p+2)])
       colnames(mat1)<-colnames(mat2)<-c("Estimates","StdErr","Wald")
       rownames(mat1)<-listanomi[1:(p+1)]; rownames(mat2)<-listanomi[(2+p):(q+p+2)]
       cat("Uncertainty                                           ", "\n")
       print(mat1,digits=digits)
       
       cat("=======================================================================","\n")
       
       cat("Feeling                                           ", "\n")
       print(mat2,digits=digits)
       cat("=======================================================================","\n")
       
       cat("Overdispersion                                           ", "\n")
       
       mat3<-cbind(stime[(p+q+3):(p+q+s+3)],StdErr[(p+q+3):(p+q+s+3)],Wald[(p+q+3):(p+q+s+3)])
       colnames(mat3)<-c("Estimates","StdErr","Wald")
       rownames(mat3)<-listanomi[(p+q+3):(p+q+s+3)]
       print(mat3,digits=digits)
       
       x$uncertainty<-mat1
       x$feeling<-mat2
       x$overdispersion<-mat3
     }
     
     
   }
   
   
   if (family == "IHG" | family =="CUSH"){
     
     if (length(stime)>1){
         flagcov<-1
     }
     
     matout<-cbind(stime,StdErr,Wald)
     
     colnames(matout)<-c("Estimates","StdErr","Wald")
     rownames(matout)<-listanomi
     
     print(matout,digits=digits)
     
     if (family=="IHG"){
       x$preference<-matout
     } else {
       x$shelter<-matout
     }
     
   }
   
   
   # for(i in 1:np){
   #   cat(nomi[i],"      ",stime[i],"        ",errstd[i],"  ",wald[i],"      ","\n")
   # }
   cat("=======================================================================","\n")
   if (family=="CUB" & !is.null(object$ellipsis$shelter)){
     covshe<-model.matrix(modello,data=mf,rhs=3)
     
     if (ncol(covshe)==0){
       X<-NULL
     } else {
       X<-covshe[,-1]
     }

     
     if (is.null(X)){
       matout<-cbind(stime,StdErr,Wald)
       
       colnames(matout)<-c("Estimates","StdErr","Wald")
       rownames(matout)<-listanomi
       print(matout,digits=digits)
       
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
       Wald<-wald2
       matout2<-cbind(mat2,ErrStd,wald2)
       dimnames(matout2)<-list(nomi2,c("Estimates","StdErr","Wald"))
       
       #matout2<-cbind(nomi2,stime2,errstd2,wald2)
       #  dimnames(matout2)<-list(rep("",length(nomi2)),c("Parameters",  "ML-estimates" , "Std. err.", "Est./Std.err (Wald test)"))
       cat("=======================================================================","\n")
       
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
   if (flagcov==0){
     cat("Log-lik(UNIFORM)   =",round(x$llunif,digits=digits),"\n")
     cat("Log-lik(saturated) =",round(x$logsat,digits=digits),"\n")
     cat("Deviance           =",round(x$devian,digits=digits),"\n")
   }
   cat("-----------------------------------------------------------------------","\n")
   #cat("Log-lik(Shifted-BINOMIAL) =",round(llsb,digits=8),"\n")
   #cat("-----------------------------------------------------------------------","\n")
   cat("AIC       =",round(x$AIC,digits=digits),"\n")
   cat("BIC       =",round(x$BIC,digits=digits),"\n")
   cat("ICOMP     =",round(x$ICOMP,digits=digits),"\n")
   cat("=======================================================================","\n")
   cat("Elapsed time=",object$time,"seconds","=====>>>",date(),"\n")
   cat("=======================================================================","\n")
   
   class(x)<-"summary.GEM"
   #return(list(x$uncertainty,x$feeling,x$overdispersion,x$shelter,x$preference))
 }  ## chiude definizione print.summary
 
 print(output)
 invisible(output$results)
  
   #invisible(output)

}

