#' @title Plot an estimated CUB model with shelter
#' @description Plotting facility for the CUB estimation of ordinal responses when a shelter effect is included
#' @aliases cubshevisual
#' @usage cubshevisual(ordinal,shelter,...)
#' @param ordinal Vector of ordinal responses (factor type)
#' @param shelter Category corresponding to the shelter choice
#' @param ... Additional arguments to be passed to plot() and text()
#' @details It represents an estimated CUB model with shelter effect as a point
#'  in the parameter space with shelter parameter indicated as label. 
#' @return A plot of the estimated parameter vector \eqn{(\pi, \xi)} as a point in the parameter space, 
#'  labeled with the shelter parameter \eqn{\delta}.
#' @keywords device
#' @export cubshevisual
#' @seealso \code{\link{cubvisual}}, \code{\link{multicub}}
#' @import graphics
#' @examples
#' data(univer)
#' attach(univer)
#' cubshevisual(global,shelter=7,digits=3,col="blue")


cubshevisual<-function(ordinal,shelter,...){
  
  
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

  if (!is.factor(ordinal)){
    stop("Response must be an ordered factor")
  }
  
 # ordinal<-factor(ordinal,ordered=TRUE)
  lev<-levels(ordinal)
  m<-length(lev)
  
  F0<-Formula(ordinal~0|0|0)
  data<-as.data.frame(ordinal)
  
  mod<-GEM(F0,data=data,shelter=shelter,family="cub")
  stime<-mod$estimates
  pai1<-stime[1];pai2<-stime[2];csi<-stime[3]
  deltaval<-round(1-pai1-pai2,digits=digits)
  paistar<-pai1/(pai1+pai2)
  plot(1-paistar,1-csi,pch=pch,col=col,main="",xlim=xlim,ylim=ylim,xlab=expression(paste("Uncertainty  ", (1-pi))),
       ylab=expression(paste("Feeling  ", (1-xi))))
  
  text(1-paistar,1-csi, labels=bquote(delta == .(deltaval)),pos=pos,offset=offset,font=font,cex=cex,col=col)
  
 
}