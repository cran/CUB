#' @title S3 method  vcov()for class "GEM"
#' @description S3 method: vcov for objects of class \code{\link{GEM}}. 
#' @aliases vcov.GEM
#' @param object An object of class \code{\link{GEM}}
#' @param ...  Other arguments
#' @method vcov GEM
#' @export 
#' @return Variance-covariance matrix of the final ML estimates for parameters of the fitted GEM model.
#' It returns the square of the estimated standard error for CUSH and IHG models with no covariates.
#' @import methods
#' @seealso \code{\link{varmatCUB}}, \code{\link{varmatCUBE}}, \code{\link{GEM}}
#' @keywords package
#' @rdname vcov.GEM



#vcov <- function(object,...) UseMethod("vcov", object)


vcov.GEM<-function(object, ...){
  
  arguments<-list(...)
  
  digits<-arguments$digits
  
  if (is.null(digits)){
    digits<-options()$digits
  }
  
  ellipsis<-object$ellipsis
  family<-object$family
  varcov<-as.matrix(object$varmat)
  listanomi<-parnames(object)
  
  if (NROW(varcov)>1){
    dimnames(varcov)<-list(listanomi,listanomi)
  } else {
    dimnames(varcov)<-list(listanomi,"Squared Standard Error")
  }
  return(round(varcov,digits=digits))
}
