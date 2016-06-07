#' @title Joint plot of estimated CUB models in the parameter space
#' @description Return a plot of estimated CUB models represented as points in the parameter space.
#' @aliases multicub
#' @usage multicub(listord, mlist, labelpoints=as.character(1:length(mlist)),
#' caption="CUB models", colours=rep("black",length(mlist)), symbols=19,
#' thickness=1.5, xwidth=c(0,1), ywidth=c(0,1), pos=rep(3,length(mlist)))
#' @export multicub
#' @param listord List of vectors of ordinal observations, possibly with different length and over different
#' numbers of categories
#' @param mlist List of numbers of categories corresponding to ordinal data in listord
#' @param labelpoints Character strings indicating the labels for the estimated models as points
#'  in the parameter space (default is the corresponding index of column in matord)
#' @param caption Character string indicating the plot title (default is "CUB models")
#' @param colours Indicate the colours for the plotted points
#' @param symbols Indicate the symbols used to represent estimated CUB models as points in the
#'  parameter space
#' @param thickness Indicate the thickness
#' @param xwidth Indicate the width of the abscissa axis
#' @param ywidth Indicate the width of the ordinate axis
#' @param pos Vector of positions relative to location for labelpoints
#' @keywords device
#' @examples
#' data(univer)
#' listord<-as.list(univer[,8:12])
#' mlist<-as.list(rep(7,length(listord)))
#' multicub(listord, mlist, labelpoints=as.character(1:length(mlist)),
#' caption="CUB models", colours=rep("black",length(mlist)), symbols=19,
#' thickness=1.5, xwidth=c(0,1), ywidth=c(0,1), pos=rep(3,length(mlist)))
#' ###############################
#' m1<-5; m2<-7;  m3<-9
#' mlist<-list(m1,m2,m3)
#' pai<-0.7;csi<-0.6
#' n1<-1000; n2<-500; n3<-1500
#' ord1<-simcub(n1,m1,pai,csi)
#' ord2<-simcub(n2,m2,pai,csi)
#' ord3<-simcub(n3,m3,pai,csi)
#' listord<-list(ord1,ord2,ord3)
#' multicub(listord,mlist,labelpoints=c("m=5","m=7","m=9"),pos=c(3,1,4))


multicub<-function(listord,mlist,labelpoints=as.character(1:length(mlist)),
                    caption="CUB models",colours=rep("black",length(mlist)), symbols=19,
                    thickness=1.5,xwidth=c(0,1),ywidth=c(0,1),pos=rep(3,length(mlist))){
  
  k<-length(listord)
  vettpai<-vettcsi<-rep(NA,k);
  plot(c(0,1),c(0,1),main=caption,cex=1.2,cex.main=1, font.lab=4,cex.lab=1,
       pch=symbols,lwd=thickness,las=1,type="n",
       xlim=xwidth,ylim=ywidth,
       xlab=expression(paste("Uncertainty  ", (1-pi))),
       ylab=expression(paste("Feeling  ", (1-xi))));
  for(j in 1:k){
    m<-as.numeric(mlist[j])
    ord<-listord[[j]]
    stimacub <- cubforsim(m,ord)
    param<-stimacub$estimates
    vettpai[j]<-param[1]; vettcsi[j]<-param[2];
    points(1-vettpai[j],1-vettcsi[j],col=colours[j],pch=19)
    
  }
  text(1-vettpai,1-vettcsi,labels=labelpoints,pos=pos,offset=0.5,font=4,cex=0.7,col=colours)
  
}
