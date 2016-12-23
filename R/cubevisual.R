#' @title Plot an estimated CUBE model
#' @description Plotting facility for the CUBE estimation of ordinal responses. 
#' @aliases cubevisual
#' @usage cubevisual(ordinal,...)
#' @param ordinal Vector of ordinal responses (factor type)
#' @param ... Additional arguments to be passed to \code{plot()} and \code{text()}
#' @details It represents an estimated CUBE model as a point
#'  in the parameter space with the overdispersion being labeled. 
#' @return A plot of the estimated parameter vector \eqn{(\pi, \xi)} as a point in the parameter space 
#' with the overdispersion \eqn{\phi} being labeled.
#' @keywords device
#' @export cubevisual
#' @import graphics
#' @examples
#' data(univer)
#' attach(univer)
#' cubevisual(global,xlim=c(0,0.5),ylim=c(0.5,1),cex=0.8,digits=3,col="red")


cubevisual<-function(ordinal,...){
  
  if (!is.factor(ordinal)){
    stop("Response must be an ordered factor")
  }
  
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
    main<-"CUBE models parameter space"
  }
    
  m<-length(levels(ordinal))
  starting<-inibestcube(m,ordinal)
  
  F0<-Formula(ordinal~0|0|0)
    stimacube<-GEM(F0,family="cube",starting=starting,maxiter = 500, toler = 1e-06)
  param<-stimacube$estimates; pai<-param[1];csi<-param[2];phi<-param[3]
  plot(1-pai,1-csi,main=main,las=1,pch=pch,cex=cex,xlim=xlim,ylim=ylim,
       col=col,
       xlab=expression(paste("Uncertainty  ", (1-pi))),
       ylab=expression(paste("Feeling  ", (1-xi))));
  text(1-pai,1-csi,labels=bquote(phi == .(round(phi,digits=3))),font=font,pos=pos,offset=offset,cex=cex,col=col)
}


