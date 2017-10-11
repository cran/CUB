#' @title Simulation routine for CUB models
#' @aliases simcub
#' @description Generate \eqn{n} pseudo-random observations following the given CUB distribution.
#' @keywords distribution
#' @usage simcub(n,m,pai,csi)
#' @export simcub
#' @import stats
#' @param n Number of simulated observations
#' @param m Number of ordinal categories
#' @param pai Uncertainty parameter
#' @param csi Feeling parameter
#' @seealso \code{\link{probcub00}}
#' @examples
#' n<-300
#' m<-9
#' pai<-0.4
#' csi<-0.7
#' simulation<-simcub(n,m,pai,csi)
#' plot(table(simulation),xlab="Ordinal categories",ylab="Frequencies")



simcub <-
function(n,m,pai,csi){
  dico<-runif(n)<pai;
  vett<-dico*(1+rbinom(n,m-1,1-csi))+(1-dico)*sample(m,n,replace=TRUE)
  return(vett)
}
