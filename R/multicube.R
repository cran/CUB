#' @title Joint plot of estimated CUBE models in the parameter space
#' @description Return a plot of estimated CUBE models represented as points in the
#'  parameter space, where the overdispersion is labeled.
#' @aliases multicube
#' @usage multicube(listord,mvett,csiplot=FALSE,paiplot=FALSE,...)
#' @export multicube
#' @param listord A data matrix, data frame, or list of vectors of ordinal observations (for variables 
#'  with different number of observations)
#' @param mvett Vector of number of categories for ordinal variables in \code{listord} (optional: if missing, 
#' the number of categories is retrieved from data: it is advisable to specify it in case some category has zero 
#' frequency)
#' @param csiplot Logical: should \eqn{\xi} or \eqn{1-\xi} be the \eqn{y} coordinate
#' @param paiplot Logical: should \eqn{\pi} or \eqn{1-\pi} be the \eqn{x} coordinate
#' @param ... Additional arguments to be passed to \code{\link{plot}}, \code{\link{text}}, and \code{\link{GEM}}
#' @keywords device
#' @return Fit a CUBE model to list elements, and then by default it returns a plot of the estimated 
#' \eqn{(1-\pi, 1-\xi)} as points in the parameter space, labeled with the estimated overdispersion.
#'  Depending on \code{csiplot} and \code{paiplot} and on desired output, \eqn{x} and \eqn{y} 
#'  coordinates may be set to \eqn{\pi} and \eqn{\xi}, respectively.
#' @examples
#' m1<-5; m2<-7;  m3<-9
#' pai<-0.7;csi<-0.6;phi=0.1
#' n1<-1000; n2<-500; n3<-1500
#' ord1<-simcube(n1,m1,pai,csi,phi)
#' ord2<-simcube(n2,m2,pai,csi,phi)
#' ord3<-simcube(n3,m3,pai,csi,phi)
#' listord<-list(ord1,ord2,ord3)
#' multicube(listord,labels=c("m=5","m=7","m=9"),pos=c(3,1,4),expinform=TRUE)


multicube<-function(listord,mvett,csiplot=FALSE,paiplot=FALSE,...){
  
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
    main<-"CUBE models"
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
  

  vettpai<-vettcsi<-vettphi<-rep(NA,k);

  labelpoints<-c()
  
  for(j in 1:k){
    
    ord<-listord[[j]]
  
    m<-length(levels(factor(ord,ordered=TRUE)))
    starting<-inibestcube(m,ord)
    F0<-Formula(ord~0|0|0)
    data<-as.data.frame(ord)
   stimacube <- GEM(F0,data=data, family="cube",starting=starting,m=mvett[j])
    #stimacube <- CUBE(F0,data=data, family="cube",starting=starting,m=listm[j])
    param<-stimacube$estimates
    vettpai[j]<-param[1]; vettcsi[j]<-param[2]; vettphi[j]<-round(param[3],digits=3)
    labelpoints[j]=as.character(paste("phi=",vettphi[j]))
    
  }
  labels<-ellipsis.arg$labels
  if(is.null(labels)){
    labels<-labelpoints
  }
  plot(c(0,1),c(0,1),main=main,cex=cex,cex.main=1, font.lab=4,cex.lab=1,
       pch=pch,las=1,type="n",
       xlim=xlim,ylim=ylim,
       xlab=xlab,ylab=ylab);
  
  
  paival<-1-vettpai; csival<-1-vettcsi
  if (csiplot==TRUE){
    csival<-vettcsi
  } 
  if (paiplot==TRUE){
    paival<-vettpai
  } 
  
  points(paival,csival,col=colours,pch=pch,cex=cex)
  text(paival,csival,labels=labels,pos=pos,offset=offset,font=font,cex=cex,
       col=colours)
 
  
}