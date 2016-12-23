#' @title  S3 method "fitted" for class "GEM"
#' @description S3 method fitted for objects of class \code{\link{GEM}}. 
#' @aliases fitted.GEM
#' @param object An object of class \code{\link{GEM}}
#' @param ...  Other arguments
#' @method fitted GEM
#' @export 
#' @details Returns the fitted probability distribution for GEM models with no covariates. If only one dichotomous 
#' covariate is included in the model to explain some components, it returns the fitted probability distribution for each profile.
#' @import methods
#' @rdname fitted.GEM
#' @keywords package
#' @examples
#' data(univer)
#' attach(univer)
#' fitcub<-GEM(Formula(global~0|freqserv|0),family="cub")
#' fitted(fitcub,digits=4)

#fitted <- function(object,...) UseMethod("fitted", object)


fitted.GEM<-function(object, ...){
  
  arguments<-list(...)
  
  digits<-arguments$digits
  
  if (is.null(digits)){
    digits<-options()$digits
  }

  theorpr<-profiles(object)
  return(round(theorpr,digits=digits))
}



profiles <- function(object) UseMethod("profiles", object)


profiles.CUB<-function(object){
  
  ellipsis<-object$ellipsis
  m<-ellipsis$m
  
  values<-c()
  for (j in 1:m){
    values[j]<-paste("R =",j)
  }
  
  stime<-object$estimates
  
  modello<-object$formula
  data<-object$data
  
  
  covpai<-model.matrix(modello,data=data,rhs=1)
  covcsi<-model.matrix(modello,data=data,rhs=2)
  covshe<-model.matrix(modello,data=data,rhs=3)
  
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
  
  
  
  if (!is.null(ellipsis$shelter)){  
    shelter<-ellipsis$shelter
    #nomi<-rbind("pai1","pai2","csi")
    pai1<-stime[1];pai2<-stime[2];csi<-stime[3]
    delta<-1-pai1-pai2
    paistar<-pai1/(pai1+pai2)
    
    nprof<-1
    theorpr<-matrix(NA,nrow=m,ncol=nprof)
    theorpr[,1]<-probcubshe1(m,pai1,pai2,csi,shelter)
    profili<-""
    
    dimnames(theorpr)<-list(values,profili)
    return(theorpr)
    
  } else {
    if ( is.null(W) & is.null(Y)){
      #  nomi<-rbind("pai","csi");
      pai<-stime[1];csi<-stime[2];
      
      nprof<-1
      theorpr<-matrix(NA,nrow=m,ncol=nprof)
      theorpr[,1]<-probcub00(m,pai,csi)
      profili<-""
      dimnames(theorpr)<-list(values,profili)
      return(theorpr)
      
    } 
    if (is.null(W) & !is.null(Y)) {
      bet<-stime[1:(length(stime)-1)]; csi<-stime[length(stime)]
      #nomi<-c(paste("beta",0:(length(bet)-1),sep="_"),"csi   ")
      #Y<-as.matrix(ellipsis$Y)
      if (NCOL(Y)==1 && length(unique(Y))==2) {
        
        Y<-as.matrix(Y)
        ny<-NCOL(Y)
        yval<-list()
        for (j in 1:ny){
          yval[[j]]<-sort(unique(Y[,j]))
        }
        
        profiles<-expand.grid(yval)
        nprof<-NROW(profiles)
        
        theorpr<-matrix(NA,nrow=m,ncol=nprof)
        
        #paivett<-unique(logis(Y,bet))
        paivett<-c()
        profili<-c()
        for (j in 1:nprof){
          profili[j]<-"("
          paivett[j]<-1/(1+ exp(-bet[1]-sum(as.numeric(profiles[j,])*bet[2:length(bet)])))
          
          theorpr[,j]<-probcub00(m,paivett[j],csi)
          for (k in 1:NCOL(profiles)){
            profili[j]<-paste(profili[j],"Y",k,"=",profiles[j,k],"")
          }
          profili[j]<-paste(profili[j],")")
        }
        dimnames(theorpr)<-list(values,profili)
        return(theorpr)
      } else {
        cat("No fitted method available","\n")
      }
      
    } 
    
    if (!is.null(W) & is.null(Y)){
    #  W<-as.matrix(ellipsis$W)
      
      if (NCOL(W)==1 && length(unique(W))==2){
        pai<-stime[1]; gama<-stime[2:length(stime)];
        # nomi<-c("pai    ",paste("gamma",0:(length(gama)-1),sep="_"))
        wval<-list()
        
        nw<-NCOL(W)
        W<-as.matrix(W)
        for (j in 1:nw){
          wval[[j]]<-sort(unique(W[,j]))
        }
        profiles<-expand.grid(wval)
        nprof<-NROW(profiles)
        theorpr<-matrix(NA,nrow=m,ncol=nprof)
        csivett<-c()
        #csivett<-unique(logis(W,gama))
        profili<-c()
        for (j in 1:nprof){
          profili[j]<-"("
          csivett[j]<-1/(1+ exp(-gama[1]-sum(as.numeric(profiles[j,])*gama[2:length(gama)])))
          theorpr[,j]<-probcub00(m,pai,csivett[j])
          for (k in 1:NCOL(profiles)){
            profili[j]<-paste(profili[j],"W",k,"=",profiles[j,k])
          }
          profili[j]<-paste(profili[j],")")
        }
        dimnames(theorpr)<-list(values,profili)
        return(theorpr)
        
      }  else {
        cat("No fitted method available","\n")
        
      }
    }
      
    if (!is.null(Y) & !is.null(W)){
     # Y<-as.matrix(ellipsis$Y)
    #  W<-as.matrix(ellipsis$W)
      ny<-NCOL(Y)
      nw<-NCOL(W)
      
      if (ny==1 & nw==1 & length(unique(Y))==2 & length(unique(W))==2 ){
        
        if (all(Y==W)){
          
          bet<-stime[1:(ny+1)];gama<-stime[(ny+2):length(stime)];
          #nomi<-c(paste("beta",0:(length(bet)-1),sep="_"),paste("gamma",0:(length(gama)-1),sep="_"))
          
          listW<-list()
          listY<-list()
          Y<-as.matrix(Y)
          W<-as.matrix(W)
          for (j in 1:ny){
            listY[[j]]<-Y[,j]
          }
          
          for (j in 1:nw){
            listW[[j]]<-W[,j]
          }
          
          eqy<-eqw<-c()
          for (j in 1:length(listY)){
            for (k in 1:length(listW)){
              if (all(listY[[j]]==listW[[k]])){
                eqy<-c(eqy,j)
                eqw<-c(eqw,k)
              }
            }
          }
          comcov<-cbind(eqy,eqw)
          
          YW<-as.matrix(unique(t(cbind(Y,W))))
          #unique restituisce le righe di un array tolte le ripetizioni
          
          paivett<-csivett<-c()
          ywval<-list()
          for (j in 1:NCOL(t(YW))){
            ywval[[j]]<-sort(unique(t(YW)[,j]))
          }
          
          
          profiles<-unique(expand.grid(ywval))
          nprof<-NROW(profiles)
          theorpr<-matrix(NA,nrow=m,ncol=nprof)
          profili<-c()
          # csivett<-unique(logis(W,gama))
          # paivett<-unique(logis(Y,bet))
          
          for (j in 1:nprof){
            vett<-rep(NA,nw)
            if (NROW(comcov)!=0){
              vett[eqw]<-as.numeric(profiles[j,eqy])
              if (nw - NROW(comcov) > 0){
                vett[-eqw]<-as.numeric(profiles[j,(ny+1):NCOL(profiles)])
              }
            } else {
              vett<-as.numeric(profiles[j,(ny+1):NCOL(profiles)])
            }
            profili[j]="("
            for (k in 1:ny){
              profili[j]<-paste(profili[j],"Y",k,"=",profiles[j,k],", ")
            }
            for (k in 1:nw){
              profili[j]<-paste(profili[j],"W",k,"=",vett[k],", ")
            }
            profili[j]<-paste(profili[j],")")
            
            csivett[j]<- 1/(1+ exp(-gama[1]-sum(vett*gama[2:length(gama)])))
            paivett[j]<-1/(1+ exp(-bet[1]-sum(profiles[j,1:ny]*bet[2:length(bet)])))
            theorpr[,j]<-probcub00(m,paivett[j],csivett[j])
          }
          dimnames(theorpr)<-list(values,profili)
          return(theorpr)
          
        } else {
          cat("No fitted method available","\n")
        }
      } else {
        cat("No fitted method available","\n")
      }   
  } 
  
}
}


####################################

profiles.CUBE<-function(object){
  
  ellipsis<-object$ellipsis
  
  stime<-object$estimates
  m<-ellipsis$m
  
  modello<-object$formula
  # EFFE<-mod$Formula
  data<-object$data
  
  covpai<-model.matrix(modello,data=data,rhs=1)
  covcsi<-model.matrix(modello,data=data,rhs=2)
  covphi<-model.matrix(modello,data=data,rhs=3)
  
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
  
  
#  Y<-ellipsis$Y
 # W<-ellipsis$W
  #Z<-ellipsis$Z
  
  values<-c()
  for (j in 1:m){
    values[j]<-paste("R =",j)
  }
  
  if (is.null(Y) & is.null(W) & is.null(Z)){
    pai<-stime[1]; csi<-stime[2]; phi<-stime[3]
    
    nprof<-1
    theorpr<-matrix(NA,nrow=m,ncol=nprof)
    theorpr[,1]<- probcube(m,pai,csi,phi)
    profili<-""
    dimnames(theorpr)<-list(values,profili)
    return(theorpr)
    
  } else if  (!is.null(W) & is.null(Y) & is.null(Z)){
    
    if (NCOL(W)==1 && length(unique(W))==2){
      pai<-stime[1]; gama<-stime[2:(length(stime)-1)]; phi<-stime[length(stime)];
      # nomi<-c("pai    ",paste("gamma",0:(length(gama)-1),sep="_"))
      wval<-list()
      W<-as.matrix(W)
      nw<-NCOL(W)
      for (j in 1:nw){
        wval[[j]]<-sort(unique(W[,j]))
      }
      profiles<-expand.grid(wval)
      nprof<-NROW(profiles)
      theorpr<-matrix(NA,nrow=m,ncol=nprof)
      csivett<-c()
      #csivett<-unique(logis(W,gama))
      profili<-c()
      for (j in 1:nprof){
        profili[j]<-"("
        csivett[j]<-1/(1+ exp(-gama[1]-sum(as.numeric(profiles[j,])*gama[2:length(gama)])))
        theorpr[,j]<-probcube(m,pai,csivett[j],phi)
        for (k in 1:NCOL(profiles)){
          profili[j]<-paste(profili[j],"W",k,"=",profiles[j,k])
        }
        profili[j]<-paste(profili[j],")")
      }
      dimnames(theorpr)<-list(values,profili)
      return(theorpr)
      
    } else {
      cat("No fitted method available","\n")
    }
    
    
  } else if (!is.null(Y) & !is.null(W) & !is.null(Z)){
    cat("No fitted method available","\n")
    
  }
  
}


####################################

profiles.IHG<-function(object){
  
  ellipsis<-object$ellipsis
  m<-ellipsis$m
  
  values<-c()
  for (j in 1:m){
    values[j]<-paste("R =",j)
  }
  
  modello<-object$formula
  #EFFE<-mod$Formula
  data<-object$data
  
  covtheta<-model.matrix(modello,data=data,rhs=1)
  
  if (ncol(covtheta)==0){
    U<-NULL
  } else {
    U<-covtheta[,-1]
  }
  
  
  
  stime<-object$estimates
  
  #U<-ellipsis$U
  if (is.null(U)){
    nprof<-1
    theta<-stime[1]
    theorpr<-matrix(NA,nrow=m,ncol=nprof)
    theorpr[,1]<- probihg(m,theta)
    profili<-""
    dimnames(theorpr)<-list(values,profili)
    return(theorpr)
    
  } else {
    if (NCOL(U)==1 && length(unique(U))==2){
      nuest<-stime
      uval<-list()
      U<-as.matrix(U)
      nu<-NCOL(U)
      for (j in 1:nu){
        uval[[j]]<-sort(unique(U[,j]))
      }
      profiles<-expand.grid(uval)
      nprof<-NROW(profiles)
      theorpr<-matrix(NA,nrow=m,ncol=nprof)
      thetavett<-c()
      #csivett<-unique(logis(W,gama))
      profili<-c()
      for (j in 1:nprof){
        profili[j]<-"("
        thetavett[j]<-1/(1+ exp(-nuest[1]-sum(as.numeric(profiles[j,])*nuest[2:length(nuest)])))
        theorpr[,j]<-probihg(m,thetavett[j])
        for (k in 1:NCOL(profiles)){
          profili[j]<-paste(profili[j],"U",k,"=",profiles[j,k])
        }
        profili[j]<-paste(profili[j],")")
      }
      dimnames(theorpr)<-list(values,profili)
      return(theorpr)
      
    } else {
      cat("No fitted method available","\n")
    }
    
  }
  
  
  
}


####################################

profiles.CUSH<-function(object,...){
  
  ellipsis<-object$ellipsis
  
  stime<-object$estimates
  m<-ellipsis$m
  
 # X<-ellipsis$X
  shelter<-ellipsis$shelter
  
  values<-c()
  for (j in 1:m){
    values[j]<-paste("R =",j)
  }
  
  modello<-object$formula
  # EFFE<-mod$Formula
  data<-object$data
  
  covshe<-model.matrix(modello,data=data,rhs=1)
 if (ncol(covshe)==0){
   X<-NULL
 } else {
   X<-covshe[,-1]
 }
 
 
 
  
  if (is.null(X)){
    nprof<-1
    delta<-stime[1]
    theorpr<-matrix(NA,nrow=m,ncol=nprof)
    theorpr[,1]<- probcush(m,delta,shelter)
    profili<-""
    dimnames(theorpr)<-list(values,profili)
    return(theorpr)
    
  } else {
    if (NCOL(X)==1 && length(unique(X)==2)){
      omega<-stime
      xval<-list()
      X<-as.matrix(X)
      nx<-NCOL(X)
      for (j in 1:nx){
        xval[[j]]<-sort(unique(X[,j]))
      }
      profiles<-expand.grid(xval)
      nprof<-NROW(profiles)
      theorpr<-matrix(NA,nrow=m,ncol=nprof)
      deltavett<-c()
      profili<-c()
      for (j in 1:nprof){
        profili[j]<-"("
        deltavett[j]<-1/(1+ exp(-omega[1]-sum(as.numeric(profiles[j,])*omega[2:length(omega)])))
        theorpr[,j]<-probcush(m,deltavett[j],shelter)
        for (k in 1:NCOL(profiles)){
          profili[j]<-paste(profili[j],"X",k,"=",profiles[j,k],",")
        }
        profili[j]<-paste(profili[j],")")
      }
      dimnames(theorpr)<-list(values,profili)
      return(theorpr)
    } else {
      cat("No fitted method available","\n")
    }
  }
  
  
}
