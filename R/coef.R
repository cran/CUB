#' @title S3 Method: coef for class "GEM"
#' @description S3 method: coef for objects of class \code{\link{GEM}}. 
#' @aliases coef.GEM
#' @method coef GEM
#' @param object An object of class \code{\link{GEM}}
#' @param ...  Other arguments
#' @import methods
#' @return ML estimates of parameters of the fitted GEM model.
#' @details Returns estimated values of coefficients of the fitted model 
#' @export 
#' @rdname coef.GEM
#' @keywords package

#coef<- function(object,...) UseMethod("coef", object)


coef.GEM<-function(object,...){
  
  arguments<-list(...)
  
  digits<-arguments$digits
  
  if (is.null(digits)){
    digits<-options()$digits
  }
  
  ellipsis<-object$ellipsis
  listanomi<-parnames(object)
  mat<-round(as.matrix(object$estimates),digits=digits)
  dimnames(mat)<-list(listanomi,"")
  
  mat
}


