#' @title Joint plot of estimated CUB models in the parameter space
#' @description Return a plot of estimated CUB models represented as points in the parameter space.
#' @aliases multicub
#' @usage multicub(listord,...)
#' @export multicub
#' @param listord List of vectors of ordinal observations (of factor class), possibly with different lengths and over different
#' numbers of categories
#' @param ... Additional arguments to be passed to \code{\link{plot}}, \code{\link{text}}, and \code{\link{GEM}}
#' @keywords device
#' @examples
#' data(univer)
#' listord<-as.list(univer[,8:12])
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


multicub<-function(listord,...){
  
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
      pos<-rep(3,length(listord))
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
    } 
    cex<-ellipsis.arg$cex
    if (is.null(cex)){
      cex<-rep(0.5,length(listord))
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
    }  
    

  k<-length(listord)
  vettpai<-vettcsi<-rep(NA,k);
  plot(c(0,1),c(0,1),main=main,cex.main=1, font.lab=4,cex.lab=1,
       pch=pch,las=1,type="n",
       xlim=xlim,ylim=ylim,
       xlab=expression(paste("Uncertainty  ", (1-pi))),
       ylab=expression(paste("Feeling  ", (1-xi))));
  labelpoints<-c()
  for(j in 1:k){
    
    ord<-listord[[j]]
   
    
    F0<-Formula(ord~0|0|0)
    data<-as.data.frame(ord)
    stimacub <- GEM(F0,data=data,family="cub",maxiter=300,toler=1e-4)
    param<-stimacub$estimates
    vettpai[j]<-param[1]; vettcsi[j]<-param[2];
    points(1-vettpai[j],1-vettcsi[j],col=colours[j],pch=pch[j],cex=cex[j])
    labelpoints[j]<-j
  }
  labels<-ellipsis.arg$labels
  if(is.null(labels)){
    labels<-labelpoints
  }
  
  text(1-vettpai,1-vettcsi,labels=labels,pos=pos,offset=offset,font=font,cex=cex,col=colours)

}



