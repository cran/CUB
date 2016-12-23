#' @title Generic function for coefficient names 
#' @description Generic function for names of parameter estimates of object of class "GEM". 
#' @aliases parnames.GEM
#' @param object An object of class "GEM"
#' @import methods
#' @keywords internal
#' @seealso \code{\link{summary}}
#' @examples
#' data(univer);attach(univer)
#' model<-GEM(Formula(officeho~0|0|0),family="cub",shelter=7)
#' model 

parnames <- function(object) UseMethod("parnames", object)


parnames.GEM<-function(object){
  
  listanomi<-c()
  if (object$family=="CUB"){
    listanomi<- parnames.CUB(object)
  }
  if (object$family=="CUBE"){
    listanomi<- parnames.CUBE(object)
  }
  if (object$family=="IHG"){
    listanomi<-  parnames.IHG(object)
  }
  if (object$family=="CUSH"){
    listanomi<-  parnames.CUSH(object)
  }
  return(listanomi)
}


parnames.CUB<-function(object){
  
  effe<-object$formula
  #   EFFE<-modello$Formula
  data<-object$data
  covpai<-model.matrix(effe,data=data,rhs=1)
  covcsi<-model.matrix(effe,data=data,rhs=2)
  covshe<-model.matrix(effe,data=data,rhs=3)
  
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
  if (ncol(covshe)==0){
    X<-NULL
  } else {
    X<-covshe[,-1]
  }
  
 
  ellipsis<-object$ellipsis
  listanomi<-c()
  if (!is.null(ellipsis$shelter)){
    if (is.null(X) & is.null(Y) & is.null(W)){
      listanomi<-c("pai1","pai2","csi")
    } else if (!is.null(X) & !is.null(Y) & !is.null(W)){
      Y<-as.matrix(Y); W<-as.matrix(W); X<-as.matrix(X);
      p<-NCOL(Y);
      q<-NCOL(W); 
      s<-NCOL(X); 
      listanomi<-c(paste("beta",0:p,sep="_"),paste("gamma",0:q,sep="_"),paste("omega",0:s,sep="_"));
    }
    
  } else {
    if (is.null(Y) & is.null(W)){
      listanomi<-c("pai","csi")
    } 
     if (!is.null(Y) & is.null(W)){
      betacoef<-c()
      npar<-length(object$estimates)
      for (j in 1:(npar-1)){
        betacoef[j]<-paste("beta",j-1,sep="_")
      }
      listanomi<-c(betacoef,"csi")
    } 
    if (is.null(Y) & !is.null(W)){
      gamacoef<-c()
      npar<-length(object$estimates)
      for (j in 1:(npar-1)){
        gamacoef[j]<-paste("gamma",j-1,sep="_")
      }
      listanomi<-c("pai",gamacoef)
    }
    if (!is.null(Y) & !is.null(W)) {
      betacoef<-gamacoef<-c()
      Y<-as.matrix(Y); W<-as.matrix(W)
      ny<-NCOL(Y); nw<-NCOL(W);
      for (j in 1:(ny+1)){
        betacoef[j]<-paste("beta",j-1,sep="_")
      }
      for (j in 1:(nw+1)){
        gamacoef[j]<-paste("gamma",j-1,sep="_")
      }
      listanomi<-c(betacoef,gamacoef)
    }
    
  }
  return(listanomi)
}


#####################################################

parnames.CUBE<-function(object){
  
  ellipsis<-object$ellipsis
  
  effe<-object$formula
  #   EFFE<-modello$Formula
  data<-object$data
  covpai<-model.matrix(effe,data=data,rhs=1)
  covcsi<-model.matrix(effe,data=data,rhs=2)
  covphi<-model.matrix(effe,data=data,rhs=3)
  
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
  
  # Y<-ellipsis$Y
  # W<-ellipsis$W
  # Z<-ellipsis$Z
  
  listanomi<-c()
  
  if (is.null(Y) & is.null(W) & is.null(Z)){
    listanomi<-rbind("pai","csi","phi")
  } else if (is.null(Y) & is.null(Z) & !is.null(W)){
    W<-as.matrix(W)
    listanomi<-c("pai",paste("gamma",0:NCOL(W),sep="_"),"phi")
  } else if (!is.null(Y) & !is.null(Z) & !is.null(W)){
    W<-as.matrix(W); Y<-as.matrix(Y); Z<-as.matrix(Z);
    listanomi<-c(paste("beta",0:NCOL(Y),sep="_"),
                 paste("gamma",0:NCOL(W),sep="_"),
                 paste("alpha",0:NCOL(Z),sep="_"))
  } else {
    cat("CUBE models not available for this variables specification")
    listanomi<-c()
  }
  return(listanomi)
}

#####################################################

parnames.IHG<-function(object){
  ellipsis<-object$ellipsis
 # U<-ellipsis$U
 
 effe<-object$formula
 data<-object$data
#   EFFE<-mod$Formula
#   data<-mod$data
  
  covtheta<-model.matrix(effe,data=data,rhs=1)
  
  
  if (ncol(covtheta)==0){
    U<-NULL
  } else {
    U<-covtheta[,-1]
  }
  
   
  listanomi<-c()
  if (is.null(U)){
    listanomi<-"theta"
  } else {
    U<-as.matrix(U)
    listanomi<-paste("nu",0:NCOL(U),sep="_")
  }
  return(listanomi)
}

#####################################################

parnames.CUSH<-function(object){
  ellipsis<-object$ellipsis
  
 # X<-ellipsis$X
  effe<-object$formula
 # EFFE<-mod$Formula
  data<-object$data
  
  covshe<-model.matrix(effe,data=data,rhs=1)
  
  if (ncol(covshe)==0){
    X<-NULL
  } else {
    X<-covshe[,-1]
  }
  
  listanomi<-c()
  if (is.null(X)){
    listanomi<-"delta"
  } else {
    X<-as.matrix(X)
    listanomi<-paste("omega",0:NCOL(X),sep="_")
  }
  return(listanomi)
}

