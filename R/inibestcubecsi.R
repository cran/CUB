#' @title Preliminary estimates of parameters for CUBE models with covariates only for feeling
#' @description Compute preliminary parameter estimates of a CUBE model with covariates only for feeling, given
#'  ordinal responses. These estimates are set as initial values to start the E-M algorithm for such model.
#' @aliases inibestcubecsi
#' @usage inibestcubecsi(m,ordinal,W,starting,maxiter,toler)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses (factor type)
#' @param W Matrix of selected covariates to explain the feeling component
#' @param starting Starting values for preliminary estimation of a CUBE without covariate
#' @param maxiter Maximum number of iterations allowed for preliminary iterations
#' @param toler Fixed error tolerance for final estimates for preliminary iterations
#' @export inibestcubecsi
#' @details Obtain preliminary estimates for the uncertainty and the overdispersion parameters by short runs of EM. As to the feeling component, it considers the
#'   nested CUB model with covariates and calls \code{\link{inibestgama}} to derive initial estimates for the coefficients
#'   of the selected covariates for feeling.
#' @return A vector \code{(pai, gamaest, phi)}, where pai is the initial estimate for the uncertainty parameter, \code{gamaest}
#'  is the vector of initial estimates for the feeling component (including an intercept term in the first entry),
#'   and \code{phi} is the initial estimate for the overdispersion parameter
#' @keywords htest utilities
#' @seealso \code{\link{inibestcube}}, \code{\link{inibestcubecov}}, \code{\link{inibestgama}}
#' @examples
#' data(relgoods)
#' attach(relgoods)
#' isnacov<-which(is.na(Gender))
#' isnaord<-which(is.na(Tv))
#' na<-union(isnacov,isnaord)
#' ordinal<-Tv[-na]; W<-Gender[-na]
#' m<-10
#' starting<-rep(0.1,3)
#' ini<-inibestcubecsi(m,ordinal,W,starting,maxiter=100,toler=1e-3)
#' nparam<-length(ini)
#' pai<-ini[1]                 # Preliminary estimates for uncertainty component
#' gamaest<-ini[2:(nparam-1)]  # Preliminary estimates for coefficients of feeling covariates
#' phi<-ini[nparam]            # Preliminary estimates for overdispersion component


inibestcubecsi <-
function(m,ordinal,W,starting,maxiter,toler){
  
  if (!is.factor(ordinal)){
    stop("Response must be an ordered factor")
  }
  
  W<-as.matrix(W)
  
  if (ncol(W)==1){
    W<-as.numeric(W)
  }
  
  gamaest<-inibestgama(m,ordinal,W)
  stimacube<-GEM(Formula(ordinal~0|0|0),family="cube",starting=starting,maxiter=maxiter,toler=toler,expinform=FALSE)
  param<-stimacube$estimates 
  elle<-length(param)
  pai<-param[1]; phi<-param[elle];
  iniest<-as.matrix(c(pai,gamaest,phi))
  dimnames(iniest)<-list(c("pai",paste("gamma",0:NCOL(W),sep="_"),"phi"),"")
  return(iniest)
}
