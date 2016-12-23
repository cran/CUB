#' @title Joint plot of estimated CUBE models in the parameter space
#' @description Return a plot of estimated CUBE models represented as points in the
#'  parameter space, where the overdispersion is labeled.
#' @aliases multicube
#' @usage multicube(listord, ...)
#' @export multicube
#' @param listord List of vectors of ordinal observations, possibly with different lengths and over different
#' numbers of categories
#' @param ... Additional arguments to be passed to \code{\link{plot}}, \code{\link{text}}, and \code{\link{GEM}}
#' @keywords device
#' @examples
#' m1<-5; m2<-7;  m3<-9
#' pai<-0.7;csi<-0.6;phi=0.1
#' n1<-1000; n2<-500; n3<-1500
#' ord1<-factor(simcube(n1,m1,pai,csi,phi),ordered=TRUE)
#' ord2<-factor(simcube(n2,m2,pai,csi,phi),ordered=TRUE)
#' ord3<-factor(simcube(n3,m3,pai,csi,phi),ordered=TRUE)
#' listord<-list(ord1,ord2,ord3)
#' multicube(listord,labels=c("m=5","m=7","m=9"),pos=c(3,1,4),expinform=TRUE)


multicube<-function(listord,...){
  
  ellipsis.arg<-list(...)
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
    pos<-3
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
    pch<-19
  } 
  cex<-ellipsis.arg$cex
  if (is.null(cex)){
    cex<-0.5
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
  }
  
  
  k<-length(listord)
  vettpai<-vettcsi<-vettphi<-rep(NA,k);
  plot(c(0,1),c(0,1),main=main,cex=cex,cex.main=1, font.lab=4,cex.lab=1,
       pch=pch,lwd=1.5,las=1,type="n",
       xlim=xlim,ylim=ylim,
       xlab=expression(paste("Uncertainty  ", (1-pi))),
       ylab=expression(paste("Feeling  ", (1-xi))));
  labelpoints<-c()
  for(j in 1:k){
    ord<-listord[[j]]
    ord<-factor(ord,ordered=TRUE)
    m<-length(levels(ord))
    starting<-inibestcube(m,ord)
    F0<-Formula(ord~0|0|0)
    data<-as.data.frame(ord)
    stimacub <- GEM(F0,data=data, family="cube",starting=starting)
    param<-stimacub$estimates
    vettpai[j]<-param[1]; vettcsi[j]<-param[2]; vettphi[j]<-round(param[3],digits=3)
    points(1-vettpai[j],1-vettcsi[j],col=colours[j],pch=pch)
    
    labelpoints[j]=as.character(paste("phi=",vettphi[j]))
    
  }
  labels<-ellipsis.arg$labels
  if(is.null(labels)){
    labels<-labelpoints
  }
  
  text(1-vettpai,1-vettcsi,labels=labels,pos=pos,offset=offset,font=font,cex=cex,
       col=colours)
  
}