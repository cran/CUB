#' @title Joint plot of estimated CUB models in the parameter space
#' @description Return a plot of estimated CUB models represented as points in the parameter space.
#' @aliases multicub
#' @usage multicub(listord,mvett,csiplot=FALSE,paiplot=FALSE,...)
#' @export multicub
#' @param listord A data matrix, data frame, or list of vectors of ordinal observations (for variables 
#'  with different number of observations)
#' @param mvett Vector of number of categories for ordinal variables in \code{listord} (optional: if missing, 
#' the number of categories is retrieved from data: it is advisable to specify it in case some category has zero 
#' frequency)
#' @param csiplot Logical: should \eqn{\xi} or \eqn{1-\xi} be the \eqn{y} coordinate
#' @param paiplot Logical: should \eqn{\pi} or \eqn{1-\pi} be the \eqn{x} coordinate
#' @param ... Additional arguments to be passed to \code{\link{plot}}, \code{\link{text}}, and \code{\link{GEM}}
#' @return Fit a CUB model to list elements, and then by default it returns a plot of the estimated 
#' \eqn{(1-\pi, 1-\xi)} as points in the parameter space. Depending on \code{csiplot} and \code{paiplot} 
#' and on desired output, \eqn{x} and \eqn{y} coordinates may be set to \eqn{\pi} and \eqn{\xi}, respectively.
#' @keywords device
#' @examples
#' data(univer)
#' listord<-univer[,8:12]
#' multicub(listord,colours=rep("red",5),cex=c(0.4,0.6,0.8,1,1.2),
#'   pch=c(1,2,3,4,5),xlim=c(0,0.4),ylim=c(0.75,1),pos=c(1,3,3,3,3))
#' ###############################
#' m1<-5; m2<-7;  m3<-9
#' pai<-0.7;csi<-0.6
#' n1<-1000; n2<-500; n3<-1500
#' ord1<-simcub(n1,m1,pai,csi)
#' ord2<-simcub(n2,m2,pai,csi)
#' ord3<-simcub(n3,m3,pai,csi)
#' listord<-list(ord1,ord2,ord3)
#' multicub(listord,labels=c("m=5","m=7","m=9"),pos=c(3,1,4))


multicub<-function(listord,mvett,csiplot=FALSE,paiplot=FALSE,...){
  
    ellipsis.arg<-list(...)
    
    if (is.data.frame(listord)==TRUE){
     listord<- as.list(listord)
    } else if (is.matrix(listord)==TRUE) {
      mat<-as.data.frame(listord)
      listord<- as.list(mat)
      
    }
    
    
    k<-length(listord)
    
    if (missing(mvett)){
      mvett<-c()
      for (j in 1:k){
        
        lev <- levels(factor(listord[[j]],ordered=TRUE)); 
        mvett[j] <- length(lev) 
      }
    }
    
    xlim<-ellipsis.arg$xlim
    if (is.null(xlim)){
      xlim<-c(0,1)
    }
    ylim<-ellipsis.arg$ylim
    if (is.null(ylim)){
      ylim<-c(0,1)
    }
    pos<-ellipsis.arg$pos
    if (is.null(pos)){
      pos<-rep(3,length(listord))
    } else if (length(pos)==1){
      pos<-rep(pos,length(listord))
    }
    
    
    offset<-ellipsis.arg$offset
    if (is.null(offset)){
      offset<-0.5
    }

    font<-ellipsis.arg$font
    if(is.null(font)){
      font<-4
    }
    pch<-ellipsis.arg$pch
    if (is.null(pch)){
      pch<-rep(19,length(listord))
    } else if (length(pch)==1){
      pch<-rep(pch,length(listord))
    }
    
    cex<-ellipsis.arg$cex
    if (is.null(cex)){
      cex<-rep(0.5,length(listord))
    } else if (length(cex)==1){
      cex<-rep(cex,length(listord))
    }
    
    
    main<-ellipsis.arg$main
    if(is.null(main)){
      main<-"CUB models"
    }
    
    labels<-ellipsis.arg$labels
    if(is.null(labels)){
      labels<-as.character(1:length(listord))
    }
    colours<-ellipsis.arg$colours
    if(is.null(colours)){
      colours<-rep("black",length(listord))
    }  else if (length(colours)==1){
      colours<-rep(colours,length(listord))
    }
    
    xlab<-ellipsis.arg$xlab
    if(is.null(xlab)){
      xlab<-expression(paste("Uncertainty  ", (1-pi)))
      if (paiplot==TRUE){
        xlab<-expression(pi)
      }
    }  
    ylab<-ellipsis.arg$ylab
    if(is.null(ylab)){
      ylab<-expression(paste("Feeling  ", (1-xi)))  
      if (csiplot==TRUE){
        ylab<-expression(xi)
      }
    }

    
    
  #listm<- mget('listm',ifnotfound=list(mlist))
  #listm<-listm$listm

  vettpai<-vettcsi<-rep(NA,k);
 
  labelpoints<-c()
  for(j in 1:k){
    
    ord<-listord[[j]]
   
    
    F0<-Formula(ord~0|0|0)
    data<-as.data.frame(ord)
   #stimacub<-CUB(F0,data=data,maxiter=300,toler=1e-4,m=listm[j])
     stimacub <- GEM(F0,data=data,family="cub",maxiter=300,toler=1e-4,m=mvett[j])
    param<-stimacub$estimates
    vettpai[j]<-param[1]; vettcsi[j]<-param[2];
    labelpoints[j]<-j
  }
  
  labels<-ellipsis.arg$labels
  if(is.null(labels)){
    labels<-labelpoints
  }
  
  plot(c(0,1),c(0,1),main=main,cex.main=1, font.lab=4,cex.lab=1,
       pch=pch,las=1,type="n", xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab);
  
  paival<-1-vettpai; csival<-1-vettcsi
  if (csiplot==TRUE){
    csival<-vettcsi
  } 
  if (paiplot==TRUE){
    paival<-vettpai
  } 
  
  

  points(paival,csival,col=colours,pch=pch,cex=cex)
  text(paival,csival,labels=labels,pos=pos,offset=offset,font=font,cex=cex,col=colours)

}



