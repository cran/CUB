#' @title Plot an estimated CUB model with shelter
#' @description Plotting facility for the CUB estimation of ordinal responses when a shelter effect is included
#' @aliases cubshevisual
#' @usage cubshevisual(ordinal,shelter,csiplot=FALSE,paiplot=FALSE,...)
#' @param ordinal Vector of ordinal responses 
#' @param shelter Category corresponding to the shelter choice
#' @param csiplot Logical: should \eqn{\xi} or \eqn{1-\xi} be the \eqn{y} coordinate
#' @param paiplot Logical: should \eqn{\pi} or \eqn{1-\pi} be the \eqn{x} coordinate
#' @param ... Additional arguments to be passed to \code{plot()} and \code{text()}. Optionally, the number \code{m} 
#' of ordinal categories may be passed: this is recommended if some category has zero frequency.
#' @details It represents an estimated CUB model with shelter effect as a point
#'  in the parameter space with shelter estimate indicated as label. 
#' @return For a CUB model with shelter fitted to \code{ordinal}, by default it returns a plot of the estimated 
#' \eqn{(1-\pi, 1-\xi)} as a point in the parameter space, labeled with the estimated shelter parameter \eqn{\delta}.
#' Depending on \code{csiplot} and \code{paiplot} and on desired output, \eqn{x} and \eqn{y} coordinates may be set
#' to \eqn{\pi} and \eqn{\xi}, respectively.
#' @keywords device
#' @export cubshevisual
#' @seealso \code{\link{cubvisual}}, \code{\link{multicub}}
#' @import graphics
#' @examples
#' data(univer)
#' ordinal<-univer$global
#' cubshevisual(ordinal,shelter=7,digits=3,col="blue",main="Global Satisfaction")


cubshevisual<-function(ordinal,shelter,csiplot=FALSE,paiplot=FALSE,...){
  
  
  ellipsis.arg<-list(...)
  
  digits<-ellipsis.arg$digits
  
  if (is.null(digits)){
    digits<-options()$digits
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
  
  main<-ellipsis.arg$main
  if (is.null(main)){
    main<-"CUB-she parameter space"
  }
  
  xlim<-ellipsis.arg$xlim
  if (is.null(xlim)){
    xlim=c(0,1)
  }
  ylim<-ellipsis.arg$ylim
  if (is.null(ylim)){
    ylim<-c(0,1)
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
  
  
 
  F0<-Formula(ordinal~0|0|0)
  #data<-as.data.frame(ordinal)
  
  mod<-GEM(F0,shelter=shelter,family="cub",m=m)
  stime<-mod$estimates
  pai1<-stime[1];pai2<-stime[2];csi<-stime[3]
  deltaval<-round(1-pai1-pai2,digits=digits)
  paistar<-pai1/(pai1+pai2)
  
  valcsi<-1-csi; valpai<-1-paistar;
  
  if (csiplot==TRUE){
    valcsi<-csi
    ylab<-expression(xi)
  }
  if (paiplot==TRUE){
    valpai<-paistar
    xlab<-expression(pi)
  }
  
  
  plot(valpai,valcsi,pch=pch,col=col,main=main,xlim=xlim,ylim=ylim,xlab=xlab,
       ylab=ylab)
  
  text(valpai,valcsi,labels=bquote(delta == .(deltaval)),pos=pos,offset=offset,font=font,cex=cex,col=col)
  
 
}