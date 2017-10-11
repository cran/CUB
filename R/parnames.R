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
  data<-object$ellipsis$data
  mf<-model.frame(effe,data=data,na.action=na.omit)
  
  covpai<-model.matrix(effe,data=mf,rhs=1)
  covcsi<-model.matrix(effe,data=mf,rhs=2)
  covshe<-model.matrix(effe,data=mf,rhs=3)
  
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
    if (NCOL(covshe)==2){
      X<-as.matrix(covshe[,-1])
      colnames(X)<-colnames(covshe)[2]
    } else {
      X<-covshe[,-1]
    }
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
      
      if (is.null(colnames(Y))){
        nomiY<- paste("beta",0:p,sep="_")
      } else {
        nomiY<-c("constant",colnames(Y))
      }
      if (is.null(colnames(W))){
        nomiW<- paste("gamma",0:q,sep="_")
      } else {
        nomiW<-c("constant",colnames(W))
      }
      if (is.null(colnames(X))){
        nomiX<- paste("omega",0:s,sep="_")
      } else {
        nomiX<-c("constant",colnames(X))
      }
      
      listanomi<-c(nomiY,nomiW,nomiX);
    }
    
  } else {
    if (is.null(Y) & is.null(W)){
      listanomi<-c("pai","csi")
    } 
     if (!is.null(Y) & is.null(W)){
      betacoef<-c()
      npar<-length(object$estimates)
      if (!is.null(colnames(Y))){
        betacoef<-c("constant",colnames(Y))
      } else {
      
       for (j in 1:(npar-1)){
        betacoef[j]<-paste("beta",j-1,sep="_")
       }
      }
       listanomi<-c(betacoef,"csi")
     } 
    
    if (is.null(Y) & !is.null(W)){
      gamacoef<-c()
      npar<-length(object$estimates)
     
      if (!is.null(colnames(W))){
        gamacoef<-c("constant",colnames(W))
      } else {
        for (j in 1:(npar-1)){
          gamacoef[j]<-paste("gamma",j-1,sep="_")
        }
      }
      
      listanomi<-c("pai",gamacoef)
    }
    
    if (!is.null(Y) & !is.null(W)) {
      betacoef<-gamacoef<-c()
      Y<-as.matrix(Y); W<-as.matrix(W)
      ny<-NCOL(Y); nw<-NCOL(W);
      
      if (is.null(colnames(Y))){
        for (j in 1:(ny+1)){
          betacoef[j]<-paste("beta",j-1,sep="_")
        }
        
        } else {
          betacoef<-c("constant",colnames(Y))
        }
        
      
      
      if (is.null(colnames(W))){
        for (j in 1:(nw+1)){
          gamacoef[j]<-paste("gamma",j-1,sep="_")
        }
        
      } else {
        gamacoef<-c("constant",colnames(W))
        
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
  data<-object$ellipsis$data
  mf<-model.frame(effe,data=data,na.action=na.omit)
  
  covpai<-model.matrix(effe,data=mf,rhs=1)
  covcsi<-model.matrix(effe,data=mf,rhs=2)
  covphi<-model.matrix(effe,data=mf,rhs=3)
  
  
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

  if (ncol(covphi)==0){
    Z<-NULL
  } else {
    if (NCOL(covphi)==2){
      Z<-as.matrix(covphi[,-1])
      colnames(Z)<-colnames(covphi)[2]
    } else {
      Z<-covphi[,-1]
    }
  }
    
  listanomi<-c()
  
  if (is.null(Y) & is.null(W) & is.null(Z)){
    listanomi<-rbind("pai","csi","phi")
  } else if (is.null(Y) & is.null(Z) & !is.null(W)){
    W<-as.matrix(W)
    
    if (is.null(colnames(W))){
      gamacoef<- paste("gamma",0:NCOL(W),sep="_")
    } else {
      gamacoef<-c("constant",colnames(W))
    }
    
   listanomi<-c("pai",gamacoef,"phi")
   
  } else if (!is.null(Y) & !is.null(Z) & !is.null(W)){
    W<-as.matrix(W); Y<-as.matrix(Y); Z<-as.matrix(Z);
    
    
    if (is.null(colnames(W))){
      gamacoef<- paste("gamma",0:NCOL(W),sep="_")
    } else {
      gamacoef<-c("constant",colnames(W))
    }
    
    if (is.null(colnames(Y))){
      betacoef<- paste("beta",0:NCOL(Y),sep="_")
    } else {
      betacoef<-c("constant",colnames(Y))
    }
    
    if (is.null(colnames(Z))){
      alfacoef<- paste("alpha",0:NCOL(Z),sep="_")
    } else {
      alfacoef<-c("constant",colnames(Z))
    }
    
    listanomi<-c(betacoef, gamacoef, alfacoef)
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
 data<-object$ellipsis$data
 mf<-model.frame(effe,data=data,na.action=na.omit)
 

  covtheta<-model.matrix(effe,data=mf,rhs=1)
  
  
  if (ncol(covtheta)==0){
    U<-NULL
  } else {
    if (NCOL(covtheta)==2){
      U<-as.matrix(covtheta[,-1])
      colnames(U)<-colnames(covtheta)[2]
    } else {
      U<-covtheta[,-1]
    }
  }
   
   
  listanomi<-c()
  if (is.null(U)){
    listanomi<-"theta"
  } else {
    U<-as.matrix(U)
    
    if (is.null(colnames(U))){
      listanomi<-paste("nu",0:NCOL(U),sep="_")
    } else {
      listanomi<-c("constant",colnames(U))
    }
  }
  return(listanomi)
}

#####################################################

parnames.CUSH<-function(object){
  ellipsis<-object$ellipsis
  
 # X<-ellipsis$X
  effe<-object$formula
 # EFFE<-mod$Formula
  data<-object$ellipsis$data
  mf<-model.frame(effe,data=data,na.action=na.omit)
 
  covshe<-model.matrix(effe,data=mf,rhs=1)
  
 if (ncol(covshe)==0){
   X<-NULL
 } else {
   if (NCOL(covshe)==2){
     X<-as.matrix(covshe[,-1])
     colnames(X)<-colnames(covshe)[2]
   } else {
     X<-covshe[,-1]
   }
 }
  
  listanomi<-c()
  if (is.null(X)){
    listanomi<-"delta"
  } else {
    X<-as.matrix(X)
    if (is.null(colnames(X))){
      listanomi<-paste("omega",0:NCOL(X),sep="_")
    } else {
      listanomi<-c("constant",colnames(X))
    }
  }
  return(listanomi)
}

