#' @title Plot an estimated CUBE model
#' @description Plotting facility for the CUBE estimation of ordinal responses. 
#' @aliases cubevisual
#' @usage cubevisual(ordinal,csiplot=FALSE,paiplot=FALSE,...)
#' @param ordinal Vector of ordinal responses 
#' @param csiplot Logical: should \eqn{\xi} or \eqn{1-\xi} be the \eqn{y} coordinate
#' @param paiplot Logical: should \eqn{\pi} or \eqn{1-\pi} be the \eqn{x} coordinate
#' @param ... Additional arguments to be passed to \code{plot()} and \code{text()}. Optionally, the number \code{m} 
#' of ordinal categories may be passed: this is recommended if some category has zero frequency.
#' @details It represents an estimated CUBE model as a point
#'  in the parameter space with the overdispersion being labeled. 
#' @return For a CUBE model fitted to \code{ordinal}, by default it returns a plot of the estimated 
#' \eqn{(1-\pi, 1-\xi)} as a point in the parameter space, labeled with the estimated overdispersion \eqn{\phi}.
#' Depending on \code{csiplot} and \code{paiplot} and on desired output, \eqn{x} and \eqn{y} coordinates may be set
#' to \eqn{\pi} and \eqn{\xi}, respectively. 
#' @keywords device
#' @export cubevisual
#' @import graphics
#' @examples
#' data(univer)
#' ordinal<-univer$global
#' cubevisual(ordinal,xlim=c(0,0.5),main="Global Satisfaction",
#'    ylim=c(0.5,1),cex=0.8,digits=3,col="red")


cubevisual<-function(ordinal,csiplot=FALSE,paiplot=FALSE,...){
  
  
  ellipsis.arg<-list(...)
  
  digits<-ellipsis.arg$digits
  
  if (is.null(digits)){
    digits<-options()$digits
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
  
  main<-ellipsis.arg$main
  if(is.null(main)){
    main<-"CUBE parameter space"
  }
  xlab<-ellipsis.arg$xlab
  if (is.null(xlab)){
    xlab<-expression(paste("Uncertainty  ", (1-pi)))
  } 
  ylab<-ellipsis.arg$ylab
  if (is.null(ylab)){
    ylab<-expression(paste("Feeling  ", (1-xi)))
  } 
  
  m<-ellipsis.arg[['m']]
  if (is.null(m)){
    ord<-factor(ordinal,ordered=TRUE)
    lev<-levels(ord)
    m<-length(lev)
    
  }
  starting<-inibestcube(m,ordinal)
  
  F0<-Formula(ordinal~0|0|0)
  stimacube<-GEM(F0,family="cube",m=m,starting=starting,maxiter = 500, toler = 1e-06)
  
  param<-stimacube$estimates; pai<-param[1];csi<-param[2];phi<-param[3]
  
  
  valcsi<-1-csi; valpai<-1-pai;
  
  if (csiplot==TRUE){
    valcsi<-csi
    ylab<-expression(xi)
  }
  if (paiplot==TRUE){
    valpai<-pai
    xlab<-expression(pi)
  }
  
  plot(valpai,valcsi,main=main,las=1,pch=pch,cex=cex,xlim=xlim,ylim=ylim,
       col=col,xlab=xlab,ylab=ylab);
  text(valpai,valcsi,labels=bquote(phi == .(round(phi,digits=3))),font=font,pos=pos,offset=offset,cex=cex,col=col)
}


