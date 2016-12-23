#' @title Plot an estimated CUB model
#' @description Plotting facility for the CUB estimation of ordinal responses. 
#' @aliases cubvisual
#' @usage cubvisual(ordinal,...)
#' @param ordinal Vector of ordinal responses
#' @param ... Additional arguments to be passed to plot() and text()
#' @details It represents an estimated CUB model as a point
#'  in the parameter space with some useful options. 
#' @return A plot of the estimated parameter vector \eqn{(\pi, \xi)} as a point in the parameter space.
#' @keywords device
#' @export cubvisual
#' @import graphics
#' @examples
#' data(univer)
#' attach(univer)
#' cubvisual(global,xlim=c(0,0.5),ylim=c(0.5,1),cex=0.8)

cubvisual<-function(ordinal,...){
  
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
  col<-ellipsis.arg$col
  if (is.null(col)){
    col<-"black"
  } 
  if (!is.factor(ordinal)){
    stop("Response must be an ordered factor")
  }
  #ordinal<-factor(ordinal,ordered=TRUE)
  F0<-Formula(ordinal~0|0|0)
  data<-as.data.frame(ordinal)
  
  stimacub<-GEM(F0,data=data, maxiter = 500, toler = 1e-06,family="cub")
  param<-stimacub$estimates; pai<-param[1];csi<-param[2];
  plot(1-pai,1-csi,main="CUB models parameter space",las=1,pch=pch,cex=cex,xlim=xlim,ylim=ylim,
       col=col,
       xlab=expression(paste("Uncertainty  ", (1-pi))),
       ylab=expression(paste("Feeling  ", (1-xi))))
  text(1-pai,1-csi,labels="estim",font=font,pos=pos,offset=offset,cex=cex,col=col)
}



