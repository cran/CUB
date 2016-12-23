#' @title S3 BIC method for class "GEM"
#' @description S3 BIC method for objects of class \code{\link{GEM}}. 
#' @aliases BIC.GEM
#' @method BIC GEM
#' @param object An object of class "GEM"
#' @param ...  Other arguments
#' @export 
#' @return BIC index for the fitted model.
#' @rdname BIC.GEM
#' @import methods
#' @seealso \code{\link{logLik}}, \code{\link{GEM}}
#' @keywords package

#BIC<- function(object,...) UseMethod("BIC", object)

BIC.GEM<-function(object,...){
  
  arguments<-list(...)
  
  digits<-arguments$digits
  
  if (is.null(digits)){
    digits<-options()$digits
  }
  bic<-object$BIC
  round(bic,digits=digits)
}

