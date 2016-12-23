#' @title logLik S3 Method for class "GEM"
#' @description S3 method: logLik() for objects of class "GEM". 
#' @aliases logLik.GEM
#' @method logLik GEM
#' @param object An object of class "GEM"
#' @param ...  Other arguments
#' @export 
#' @return Log-likelihood at the final ML estimates for parameters of the fitted GEM model.
#' @import methods
#' @seealso \code{\link{loglikCUB}}, \code{\link{loglikCUBE}}, \code{\link{GEM}},  \code{\link{loglikIHG}},
#' \code{\link{loglikCUSH}}, \code{\link{BIC}}
#' @keywords package
#' @rdname logLik.GEM



#logLik <- function(object,...) UseMethod("logLik", object)

logLik.GEM<-function(object,...){
  
  arguments<-list(...)
  
  digits<-arguments$digits
  
  if (is.null(digits)){
    digits<-options()$digits
  }
  
  return(round(object$loglik,digits=digits))    
  
}