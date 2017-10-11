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
#' @seealso \code{\link{GEM}}, \code{\link{summary}}
#coef<- function(object,...) UseMethod("coef", object)


coef.GEM<-function(object,...){
  
  arguments<-list(...)
  
  digits<-arguments$digits
  output<-list()
  
  if (is.null(digits)){
    digits<-options()$digits
  }
  
  ellipsis<-object$ellipsis
  listanomi<-parnames(object)
  mat<-round(as.matrix(object$estimates),digits=digits)
  dimnames(mat)<-list(listanomi,"")
  mat
  
}

#   
# ################### prova 
# 
# 
# print.coef.GEM <- function(x,...){
# 
#   object<-x$object
#   family<-object$family
#   mat<-x$values
#   stime<-object$estimates
#   
#   listanomi<-parnames(object)
#   
#   modello<-object$formula
#   data<-object$data
#   
#   mf<-model.frame(modello,data=data,na.action=na.omit)
#   
#  
#   data<-object$data
# 
#   if ( family == "CUB"){
#     covpai<-model.matrix(modello,data=mf,rhs=1)
#     covcsi<-model.matrix(modello,data=mf,rhs=2)
#     covshe<-model.matrix(modello,data=mf,rhs=3)
#     
#     if (ncol(covpai)==0){
#       Y<-NULL
#     } else {
#       
#       if (NCOL(covpai)==2){
#         Y<-as.matrix(covpai[,-1])
#         colnames(Y)<-colnames(covpai)[2]
#       } else {
#         Y<-covcsi[,-1]
#       }
#     }
#     if (ncol(covcsi)==0){
#       W<-NULL
#     } else {
#       if (NCOL(covcsi)==2){
#         W<-as.matrix(covcsi[,-1])
#         colnames(W)<-colnames(covcsi)[2]
#       } else {
#         W<-covcsi[,-1]
#       }
#     }
#     
#     if (ncol(covshe)==0){
#       X<-NULL
#     } else {
#       X<-covshe[,-1]
#     }
#     
#     if (!is.null(X) & !is.null(Y) & !is.null(W) & !is.null(object$ellipsis$shelter)){
#       Y<-as.matrix(Y); W<-as.matrix(W); X<-as.matrix(X);
#       p<-NCOL(Y);
#       q<-NCOL(W); 
#       s<-NCOL(X); 
#       
#       mat1<-cbind(mat[1:(p+1),1])
#       colnames(mat1)<-c("Estimates")
#       rownames(mat1)<-listanomi[1:(p+1)]
#       
#       mat2<-cbind(mat[(p+2):(p+q+2),1])
#       rownames(mat2)<-listanomi[(p+2):(p+q+2)]
#       
#       mat3<-cbind(mat[(p+q+3):(p+q+s+3),1])
#       rownames(mat3)<-listanomi[(p+q+3):(p+q+s+3)]
#       
#       cat("Uncertainty                                           ", "\n")
#       
#       print(mat1,digits=digits)
#       cat("==================================","\n")
#       cat("Feeling                                           ", "\n")
#       colnames(mat2)<-c("Estimates")
#       
#       print(mat2,digits=digits)
#       cat("==================================","\n")
#       
#       cat("Shelter effect                                ", "\n")
#       colnames(mat3)<-"Estimates"
#       print(mat3,digits=digits)
#       
#     } else if (is.null(object$ellipsis$shelter) & is.null(X) & !is.null(Y) & !is.null(W)){
#       Y<-as.matrix(Y); W<-as.matrix(W);
#       p<-NCOL(Y);
#       q<-NCOL(W); 
#       
#       mat1<-as.matrix(mat[1:(p+1),1])
#       colnames(mat1)<-c("Estimates")
#       rownames(mat1)<-listanomi[1:(p+1)]
#       
#       mat2<-as.matrix(mat[(p+2):(p+q+2),1])
#       rownames(mat2)<-listanomi[(p+2):(p+q+2)]
#       
#       cat("Uncertainty                                           ", "\n")
#       
#       print(mat1,digits=digits)
#       colnames(mat2)<-c("Estimates")
#       cat("==================================","\n")
#       
#       cat("Feeling                                           ", "\n")
#       print(mat2,digits=digits)
#       
#     } else if (is.null(object$ellipsis$shelter) & is.null(X) & is.null(Y) & !is.null(W)){
#       W<-as.matrix(W);
#       
#       q<-NCOL(W); 
#       
#       mat1<-as.matrix(mat[1,1])
#       colnames(mat1)<-c("Estimates")
#       rownames(mat1)<-listanomi[1]
#       
#       mat2<-cbind(mat[2:(q+2),1])
#       colnames(mat2)<-c("Estimates")
#       rownames(mat2)<-listanomi[2:(q+2)]
#       
#       cat("Uncertainty                                           ", "\n")
#       
#       print(mat1,digits=digits)
#       cat("==================================","\n")
#       
#       cat("Feeling                                           ", "\n")
#       print(mat2,digits=digits)
#       
#     } else if (is.null(object$ellipsis$shelter) & is.null(X) & !is.null(Y) & is.null(W)){
#       Y<-as.matrix(Y);
#       
#       p<-NCOL(Y); 
#       
#       mat1<-cbind(mat[1:(p+1),1])
#       colnames(mat1)<-c("Estimates")
#       rownames(mat1)<-listanomi[1:(p+1)]
#       
#       mat2<-cbind(mat[(p+2),1])
#       rownames(mat2)<-listanomi[p+2]
#       
#       colnames(mat2)<-c("Estimates")
#       
#       cat("Uncertainty                                           ", "\n")
#       
#       print(mat1,digits=digits)
#       cat("==================================","\n")
#       
#       cat("Feeling                                           ", "\n")
#       print(mat2,digits=digits)
#       
#     } else if (is.null(object$ellipsis$shelter) & is.null(X) & is.null(Y) & is.null(W)) {
#       
#       mat1<-cbind(mat[1])
#       colnames(mat1)<-c("Estimates")
#       rownames(mat1)<-listanomi[1]
#         
#         
#       mat2<-cbind(mat[2])
#       colnames(mat2)<-c("Estimates")
#       rownames(mat2)<-listanomi[2]
#       
#       cat("Uncertainty                                           ", "\n")
#       
#       print(mat1,digits=digits)
#       cat("==================================","\n")
#       
#       cat("Feeling                                           ", "\n")
#       print(mat2,digits=digits)
#       
#       
#     }
#   }
#   
#   
#   if (family == "CUBE"){
#     
#     covpai<-model.matrix(modello,data=mf,rhs=1)
#     covcsi<-model.matrix(modello,data=mf,rhs=2)
#     covphi<-model.matrix(modello,data=mf,rhs=3)
#     
#     if (ncol(covpai)==0){
#       Y<-NULL
#     } else {
#       Y<-covpai[,-1]
#     }
#     if (ncol(covcsi)==0){
#       W<-NULL
#     } else {
#       W<-covcsi[,-1]
#     }
#     if (ncol(covphi)==0){
#       Z<-NULL
#     } else {
#       Z<-covphi[,-1]
#     }
#     
#     if (is.null(Y)& is.null(W) & is.null(Z)){
#       mat1<-cbind(mat[1])
#       colnames(mat1)<-c("Estimates")
#       rownames(mat1)<-listanomi[1]
#       
#       mat2<-cbind(mat[2])
#       colnames(mat2)<-c("Estimates")
#       rownames(mat2)<-listanomi[2]
#       
#       cat("Uncertainty                                           ", "\n")
#       
#       print(mat1,digits=digits)
#       cat("==================================","\n")
#       
#       cat("Feeling                                           ", "\n")
#       print(mat2,digits=digits)
#       cat("==================================","\n")
#       
#       cat("Overdispersion                                           ", "\n")
#       
#       mat3<-cbind(mat[3])
#       colnames(mat3)<-c("Estimates")
#       rownames(mat3)<-listanomi[3]
#       
#       print(mat3,digits=digits)
#     } else if (is.null(Y)& !is.null(W) & is.null(Z)){
#       
#       q<-NCOL(W)
#       
#       mat1<-cbind(mat[1])
#       colnames(mat1)<-c("Estimates")
#       rownames(mat1)<-listanomi[1]
#       
#       mat2<-cbind(mat[2:(q+2)])
#       colnames(mat2)<-c("Estimates")
#       rownames(mat2)<-listanomi[2:(q+2)]
#       
#       cat("Uncertainty                                           ", "\n")
#       print(mat1,digits=digits)
#       
#       cat("==================================","\n")
#       
#       cat("Feeling                                           ", "\n")
#       print(mat2,digits=digits)
#       cat("==================================","\n")
#       
#       cat("Overdispersion                                           ", "\n")
#       
#       mat3<-cbind(mat[q+3])
#       colnames(mat3)<-c("Estimates")
#       rownames(mat3)<-listanomi[q+3]
#       
#       print(mat3,digits=digits)
#       
#       
#     }
#     else if (!is.null(Y)& !is.null(W) & !is.null(Z)){
#       p<-NCOL(Y)
#       q<-NCOL(W)
#       s<-NCOL(Z)
#       
#       mat1<-cbind(mat[1:(p+1)])
#       mat2<-cbind(mat[(2+p):(q+p+2)])
#       colnames(mat1)<-colnames(mat2)<-c("Estimates")
#       rownames(mat1)<-listanomi[1:(p+1)]
#       rownames(mat2)<-listanomi[(2+p):(q+p+2)]
#       
#       
#       cat("Uncertainty                                           ", "\n")
#       print(mat1,digits=digits)
#       
#       cat("==================================","\n")
#       
#       cat("Feeling                                           ", "\n")
#       print(mat2,digits=digits)
#       cat("==================================","\n")
#       
#       cat("Overdispersion                                           ", "\n")
#       
#       mat3<-cbind(mat[(p+q+3):(p+q+s+3)])
#       rownames(mat2)<-listanomi[(p+q+3):(p+q+s+3)]
#       
#       colnames(mat3)<-c("Estimates")
#       print(mat3,digits=digits)
#       
#       
#     }
#     
#     
#   }
#   
#   
#   if (family == "IHG" | family =="CUSH"){
#     matout<-mat
#     
#     colnames(matout)<-c("Estimates")
#     print(matout,digits=digits)
#   }
#   
#   
# 
# #  cat("==================================","\n")
#   if (family=="CUB" & !is.null(object$ellipsis$shelter)){
#     covshe<-model.matrix(modello,data=mf,rhs=3)
#     
#     if (ncol(covshe)==0){
#       X<-NULL
#     } else {
#       X<-covshe[,-1]
#     }
#     
#     
#     if (is.null(X)){
#       matout<-mat
#       
#       colnames(matout)<-c("Estimates")
#       print(matout,digits=digits)
#             pai1<-stime[1];pai2<-stime[2];csi<-stime[3]
#       delta<-1-pai1-pai2
#       paistar<-pai1/(pai1+pai2)
#       stime2<-c(paistar,csi,delta)
#       nomi2<-c("paistar","csi","delta")
#       
#       mat2<-as.matrix(stime2)
#       matout2<-mat2
#       dimnames(matout2)<-list(nomi2,c("Estimates"))
#       
#    
#       cat("==================================","\n")
#       cat("Alternative parameterization","\n")
#       print(matout2,digits=digits)
#       cat("==================================","\n")
#     }
#   }
#   
# 
# }
# print(output)





