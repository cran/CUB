#' @title Preliminary parameter estimates for CUBE models with covariates
#' @description Compute preliminary parameter estimates for a CUBE model with covariates for all the three parameters. 
#' These estimates are set as initial values to start the E-M algorithm within maximum likelihood estimation.
#' @aliases inibestcubecov
#' @usage inibestcubecov(m,ordinal,Y,W,Z)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses 
#' @param Y Matrix of selected covariates to explain the uncertainty parameter
#' @param W Matrix of selected covariates to explain the feeling parameter
#' @param Z Matrix of selected covariates to explain the overdispersion parameter
#' @export inibestcubecov
#' @return A vector \code{(inibet, inigama, inialpha)} of preliminary estimates of parameter vectors for 
#'  \eqn{\pi = \pi(\bold{\beta})}, \eqn{\xi=\xi(\bold{\gamma})}, \eqn{\phi=\phi(\bold{\alpha})}, respectively, of a CUBE model with covariates for all the three
#'   parameters. In details, \code{inibet}, \code{inigama} and \code{inialpha} have length equal to NCOL(Y)+1, NCOL(W)+1 and
#'   NCOL(Z)+1, respectively, to account for an intercept term for each component. 
#' @keywords htest utilities
#' @seealso \code{\link{inibestcube}}, \code{\link{inibestcubecsi}}, \code{\link{inibestgama}}
#' @examples
#' data(relgoods)
#' m<-10 
#' naord<-which(is.na(relgoods$Tv))
#' nacovpai<-which(is.na(relgoods$Gender))
#' nacovcsi<-which(is.na(relgoods$year.12))
#' nacovphi<-which(is.na(relgoods$EducationDegree))
#' na<-union(union(naord,nacovpai),union(nacovcsi,nacovphi))
#' ordinal<-relgoods$Tv[-na]
#' Y<-relgoods$Gender[-na]
#' W<-relgoods$year.12[-na]
#' Z<-relgoods$EducationDegree[-na]
#' ini<-inibestcubecov(m,ordinal,Y,W,Z)
#' p<-NCOL(Y)
#' q<-NCOL(W)
#' inibet<-ini[1:(p+1)]               # Preliminary estimates for uncertainty 
#' inigama<-ini[(p+2):(p+q+2)]        # Preliminary estimates for feeling 
#' inialpha<-ini[(p+q+3):length(ini)] # Preliminary estimates for overdispersion


inibestcubecov <-
function(m,ordinal,Y,W,Z){
  
  if (is.factor(ordinal)){
    ordinal<-unclass(ordinal)
  }

  Y<-as.matrix(Y); W<-as.matrix(W);Z<-as.matrix(Z)
  
  if (ncol(W)==1){
    W<-as.numeric(W)
  }
  if (ncol(Y)==1){
    Y<-as.numeric(Y)
  }
  if (ncol(Z)==1){
    Z<-as.numeric(Z)
  }
  
  q<-NCOL(Y)
  p<-NCOL(W)
  v<-NCOL(Z)
  inigama<-inibestgama(m,ordinal,W)
  inicube<-inibestcube(m,ordinal)
  pai<-inicube[1]
  bet0<-log(pai/(1-pai))
  inibet<-c(bet0,rep(0.1,q))
  alpha0<-log(0.1)
  inialpha<-c(alpha0,rep(0.1,v))
  
  ini<-as.matrix(c(inibet,inigama,inialpha))
  rownames(ini)<- c(paste("beta",0:NCOL(Y),sep="_"),
                               paste("gamma",0:NCOL(W),sep="_"),
                               paste("alpha",0:NCOL(Z),sep="_"))
  colnames(ini)<-""
  
  return(ini)
}
