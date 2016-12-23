#' @title Plot facilities for GEM objects
#' @description Plot facilities for objects of class "GEM". 
#' @aliases makeplot
#' @param object An object of class "GEM"
#' @export
#' @details Returns a plot comparing fitted 
#' probabilities and observed relative frequencies for GEM models without covariates. If only one 
#' explanatory dichotomous variable is included in the model for one or all components, 
#' then the function returns a plot comparing the distributions of the responses conditioned to
#' the value of the covariate. 
#' @keywords models device package
#' @seealso \code{\link{cubvisual}}, \code{\link{cubevisual}}, \code{\link{multicub}}, \code{\link{multicube}} 


makeplot<-function(object){
  if (object$family=="CUB"){
    makeplotCUB(object)
  }
  if (object$family=="CUBE"){
    makeplotCUBE(object)
  }
  if (object$family=="IHG"){
    makeplotIHG(object)
  }
  if (object$family=="CUSH"){
    makeplotCUSH(object)
  }
}

#makeplot <- function(object) UseMethod("makeplot", object)
# 
# makeplot.GEM<-function(object){
#   
#   if (object$family=="CUB"){
#     makeplot.CUB(object)
#   }
#   if (object$family=="CUBE"){
#     makeplot.CUBE(object)
#   }
#   if (object$family=="IHG"){
#     makeplot.IHG(object)
#   }
#   if (object$family=="CUSH"){
#     makeplot.CUSH(object)
#   }
# }

makeplotCUB<-function(object){
  
  ellipsis<-object$ellipsis
  ordinal<-object$ordinal
  family<-object$family
  m <- ellipsis$m
  
  ordinal <- unclass(ordinal)
  n<-length(ordinal)
  
  modello<-object$formula
  #EFFE<-mod$Formula
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
  
  

  stime<-round(object$estimates,5);
   
  if (!is.null(ellipsis$shelter)){
    if (is.null(W) & is.null(Y) & is.null(X)){
      theorpr<-fitted(object)
      pai1<-stime[1];pai2<-stime[2];csi<-stime[3]
      delta<-1-pai1-pai2
      paistar<-pai1/(pai1+pai2)
      
      freq<-tabulate(ordinal,nbins=m)
      dissshe<-dissim(freq/n,theorpr[,1])
      
      plot(cbind(1:m,1:m),cbind(theorpr[,1],(freq/n)),
           main=paste("CUB with shelter effect"," (Diss =",round(dissshe,digits=4),")"),
           xlim=c(1,m),ylim=c(0.0,1.1*max(theorpr[,1],(freq/n))),
           xlab="Ordinal values of R=1,2,...,m",
           ylab=expression(paste("Observed  freq. (dots) and fitted prob. (circles)")));
      
      
      points(1:m,theorpr[,1],pch=21,cex=1.2,lwd=2.0,type="b",lty=3);
      points(1:m,freq/n,pch=16,cex=1.2);
      abline(h=0);
    }  else {
      cat("No built-in plot method for this variables specifications: see multicub() and cubshevisual()","\n")
    }
  } else {
    if ( is.null(W) & is.null(Y) & is.null(X)){
      theorpr<-fitted(object)
      #freq<-matrix(NA,nrow=m,ncol=nprof)
      freq<-tabulate(ordinal,nbins=m)
      stringtitle<-"CUB model";
      thpr<-theorpr[,1]
      dissimi<-dissim(thpr,freq/n)
      #par(mfrow=c(2,1))
      plot(cbind(1:m,1:m),cbind(thpr,(freq/n)),las=1,
           main=paste(stringtitle,  "     (Diss =",round(dissimi,digits=4),")"),
           xlim=c(1,m),ylim=c(0.0,1.1*max(thpr,(freq/n))),
           xlab="Ordinal values of R=1,2,...,m",
           ylab=expression(paste("Observed relative frequencies (dots) and fitted probabilities (circles)")));
      ###
      points(1:m,thpr,pch=21,cex=1.5,lwd=2.0,type="b",lty=3); ### ex pch=8,col="red"
      points(1:m,freq/n,pch=16,cex=1.25,lwd=1.5);
      abline(h=0);
      
      # cubvisual(m,ordinal,labelpoint="estim")
      
      #  par(mfrow=c(1,1))
      
    } 
      
      if (is.null(W) & !is.null(Y) & is.null(X)) {
      
      #Y<-ellipsis$Y
      
      if (length(unique(Y))==2){
        theorpr<-fitted(object)
        prob0<-theorpr[,1]
        prob1<-theorpr[,2]
        maxpr<-max(prob0,prob1)
        
        plot(1:m,prob0,ylim=c(0.0,1.1*maxpr),cex.main=0.9,las=1,
             main="CUB distributions, given pai-covariate=0, 1",
             cex=1.2,xlab="Ordinal values of R=1,2,...,m",
             ylab="Prob(R|D=0) (circles) and  Prob(R|D=1) (dots)",pch=1,lty=1,type="b");
        lines(1:m,prob1,cex=1.2,pch=19,lty=2,type="b");
        abline(h=0);
        
      } else {
        cat("No built-in plot method for this variables specifications: see multicub() and cubvisual()","\n")
        # multicub(listaord,as.list(rep(m,nprof)),labelpoints = profili)
        
      }
      
    } 
    if (!is.null(W) & is.null(Y) & is.null(X)){
      
      #W<-as.matrix(ellipsis$W)
      
      if (NCOL(W)==1 && length(unique(W))==2){
        theorpr<-fitted(object)
        
        prob0<-theorpr[,1]
        prob1<-theorpr[,2]
        maxpr<-max(prob0,prob1)
        
        plot(1:m,prob0,ylim=c(0.0,1.1*maxpr),cex.main=0.9,las=1,
             main="CUB distributions, given csi-covariate=0, 1",
             cex=1.2,xlab="Ordinal values of R=1,2,...,m",
             ylab="Prob(R|D=0) (circles) and  Prob(R|D=1) (dots)",pch=1,lty=1,type="b");
        lines(1:m,prob1,cex=1.2,pch=19,lty=2,type="b");
        abline(h=0);
        
        #  multicub(listaord,as.list(rep(m,nprof)),labelpoints = profili)
        
      } else {
        cat("No built-in plot method for this variables specifications: see multicub() and cubvisual()","\n")
        

      }
      
    } 
    if(!is.null(W) & !is.null(Y) & is.null(X)){
      
      # Y<-as.matrix(ellipsis$Y)
      # W<-as.matrix(ellipsis$W)
      ny<-NCOL(Y)
      nw<-NCOL(W)
      
      if (ny==1 & nw==1 & length(unique(Y))==2 & length(unique(W))==2 ){
        
        if (all(Y==W)){
          theorpr<-fitted(object)
          #par(mfrow=c(2,1))
          prob0<-theorpr[,1]
          prob1<-theorpr[,2]
          maxpr<-max(prob0,prob1)
          
          plot(1:m,prob0,ylim=c(0.0,1.1*maxpr),cex.main=0.9,las=1,
               main="CUB distributions, given pai-csi covariate=0, 1",
               cex=1.2,xlab="Ordinal values of R=1,2,...,m",
               ylab="Prob(R|D=0) (circles) and  Prob(R|D=1) (dots)",pch=1,lty=1,type="b");
          lines(1:m,prob1,cex=1.2,pch=19,lty=2,type="b");
          abline(h=0);
          
          #  multicub(listaord,as.list(rep(m,nprof)),labelpoints = profili)
          #   par(mfrow=c(1,1))
        }
      } else {
        cat("No built-in plot method for this variables specifications: see multicub() and cubvisual()","\n")
      } 
      # multicub(listaord,as.list(rep(m,nprof)),labelpoints = profili)
      #par(mfrow=c(1,1))
    } 
    
  }  #chiude check su shelter
}


makeplotCUBE<-function(object){
  
  ellipsis<-object$ellipsis
  ordinal<-object$ordinal
  family<-object$family
  m <- ellipsis$m
  
  ordinal <- unclass(ordinal)
  
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
  

  
  n<-length(ordinal)
  theorpr<-fitted(object)
  
  freq<-tabulate(ordinal,nbins=m)
  
  stime<-round(object$estimates,5)
  # Y<-ellipsis$Y
  # W<-ellipsis$W
  # Z<-ellipsis$Z
  theorpr<-fitted(object)
  dissimcube<-dissim(theorpr,freq/n)
  
  if (is.null(Y) & is.null(W) & is.null(Z)){
    #par(mfrow=c(2,1))
    stringtitle="CUBE model estimation ";
    plot(cbind(1:m,1:m),cbind(theorpr,(freq/n)),las=1,
         main=paste(stringtitle,  "     (Diss =",round(dissimcube,digits=4),")"),
         xlim=c(1,m),ylim=c(0.0,1.1*max(theorpr,(freq/n))),
         xlab="Ordinal values of R=1,2,...,m",
         ylab="Obs. relative frequencies (dots) and fitted prob. (circles)",cex.lab=0.9,cex.main=0.9);
    ###
    points(1:m,theorpr,pch=21,cex=1.5,lwd=2.0,type="b",lty=3); ### ex pch=8,col="red"
    points(1:m,freq/n,pch=16,cex=1.25,lwd=1.5);
    abline(h=0);
    #
  } else if (is.null(Y) & !is.null(W) & is.null(Z)){
    
    if (NCOL(W)==1 & length(unique(W))==2){
      #par(mfrow=c(1,1))
      pai<-stime[1]; gama<-stime[2:(length(stime)-1)]; phi<-stime[length(stime)];
      vett<-as.matrix(c(0,1))
      csi0<-logis(vett[1],gama); prob0<-probcube(m,pai,csi0,phi); #theorpr[,1]? (ordine)
      csi1<-logis(vett[2],gama); prob1<-probcube(m,pai,csi1,phi);
      maxpr<-max(prob0,prob1)
      plot(1:m,prob0,ylim=c(0.0,1.1*maxpr),cex.main=0.9,las=1,
           main="CUBE distributions, given csi-covariate=0, 1",
           cex=1.2,xlab="Ordinal values of R=1,2,...,m",
           ylab="Prob(R|D=0) (circles) and  Prob(R|D=1) (dots)",pch=1,lty=1,type="b");
      lines(1:m,prob1,cex=1.2,pch=19,lty=2,type="b");
      abline(h=0);
      
    } else {
      cat("No built-in Plot method available for this variables specification: see multicube() or cubevisual()","\n")
    }
  }
}


makeplotIHG<-function(object){
  
  ellipsis<-object$ellipsis
  
  ordinal<-object$ordinal
  m <- ellipsis$m
  ordinal<-unclass(ordinal)
  freq<-tabulate(ordinal,nbins=m)
  n <-length(ordinal)
  
  modello<-object$formula
  #EFFE<-mod$Formula
  data<-object$data
  
  covtheta<-model.matrix(modello,data=data,rhs=1)
  
  
  if (ncol(covtheta)==0){
    U<-NULL
  } else {
    U<-covtheta[,-1]
  }
  
  #U<-ellipsis$U
  
  if (is.null(U)){
    theorpr<-fitted(object)
    dissihg<-dissim(theorpr[,1],freq/n)
    plot(cbind(1:m,1:m),cbind(theorpr[,1],(freq/n)),
         main=paste("IHG model (without covariates)","    (Diss =",round(dissihg,digits=4),")"),
         xlim=c(1,m),ylim=c(0,1.1*max(theorpr[,1],(freq/n))),las=1,
         xlab="Ordinal values of R=1,2,...,m",
         ylab=expression(paste("Obs. relative frequencies (dots) and fitted prob. (circles)")),
         cex.lab=0.9,cex.main=0.9)
    points(1:m,theorpr,pch=21,cex=1.5,  lwd=2.0,type="b",lty=3)
    points(1:m,freq/n,pch=16,cex=1.2)
    ### points(shelter,theorpr[shelter]-delta,pch=8);
    abline(h=0);
    #################
  } else {
    nuest<-object$estimates
    if (NCOL(U)==1 & length(unique(U))==2){
      theorpr<-fitted(object)
      vett<-as.matrix(c(0,1))
      theta0<-logis(vett[1],nuest)
      theta1<-logis(vett[2],nuest)
      prob0<-probihg(m,theta0)
      prob1<-probihg(m,theta1)
      maxpr<-max(prob0,prob1)
      
      plot(1:m,prob0,ylim=c(0.0,1.1*maxpr),cex.main=0.9,las=1,
           main="IHG distributions, given theta-covariate=0, 1",cex=1.2, cex.lab=0.9,
           xlab="Ordinal values of R=1,2,...,m",
           ylab="Prob(R|D=0) (circles) and  Prob(R|D=1) (dots)",pch=1,lty=1,type="b");
      lines(1:m,prob1,cex=1.2,pch=19,lty=2,type="b");
      abline(h=0);
      
    } else {
      cat("No built-in plot method available for this variables specification","\n")
    }
    
  }
}



makeplotCUSH<-function(object){
  
  ellipsis<-object$ellipsis
  
  ordinal<-object$ordinal
  shelter<-ellipsis$shelter
  m <- ellipsis$m
  ordinal<-unclass(ordinal)
  freq<-tabulate(ordinal,nbins=m)
  n <-length(ordinal)
  
  modello<-object$formula
 # EFFE<-mod$Formula
  data<-object$data
  
  covshe<-model.matrix(modello,data=data,rhs=1)
  
  if (ncol(covshe)==0){
    X<-NULL
  } else {
    X<-covshe[,-1]
  }
  
  #X<-ellipsis$X
  if (is.null(X)){
    theorpr<-fitted(object)
    fc<-freq[shelter]/n
    deltaest<-object$estimates
    diss00<-dissim(theorpr[,1],freq/n)
    ################### GRAFICI sovrapposti #########
    # par(mar=c(4,4,2.5,1)+0.1) ;               ### reset standard margins
    # par(mfrow=c(2,1));            ### ripristina l'area del grafico
    ###### Distributions
    stringtitle="CUSH model (without covariates)"; 
    plot(cbind(1:m,1:m),cbind(theorpr[,1],(freq/n)),las=1,
         main=paste(stringtitle," (Diss =",round(diss00,digits=4),")"),
         xlim=c(1,m),ylim=c(0.0,1.1*max(theorpr[,1],(freq/n))),
         xlab="Ordinal values of R=1,2,...,m", 
         ylab="Obs. freq (dots) and fitted prob. (circles)",
         cex.lab=0.9,cex.main=0.9); 
    points(1:m,theorpr,pch=21,cex=1.5,lwd=2.0,type="b",lty=3);
    points(1:m,freq/n,pch=16,cex=1.5,lwd=1.5);
    abline(h=0);
  } else {
    X<-as.matrix(X)
    omegaest<-object$estimates
    if (NCOL(X)==1 & length(unique(X))==2){
      theorpr<-fitted(object)
      vett<-as.matrix(c(0,1))
      delta0<-logis(vett[1],omegaest)
      delta1<-logis(vett[2],omegaest)
      prob0<-probcush(m,delta0,shelter)
      prob1<-probcush(m,delta1,shelter)
      
      maxpr<-max(prob0,prob1)
      
      plot(1:m,prob0,ylim=c(0.0,1.1*maxpr),cex.main=0.9,las=1,
           main="CUSH distributions, given delta-covariate=0, 1",cex=1.2,cex.lab=0.9,
           xlab="",ylab="Prob(R|D=0) (circles) and  Prob(R|D=1) (dots)",pch=1,lty=1,type="b");
      lines(1:m,prob1,cex=1.2,pch=19,lty=2,type="b");
      abline(h=0);
    } else {
      cat("No built-in plot method available for this variables specification","\n")
    }
    
  }
  
}



