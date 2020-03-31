#' @title Variance-covariance matrix for CUBE models
#' @aliases varmatCUBE
#' @description Compute the variance-covariance matrix of parameter estimates for CUBE models when no covariate
#' is specified, or when covariates are included for all the three parameters.
#' @usage varmatCUBE(ordinal,m,param,Y=0,W=0,Z=0,expinform=FALSE)
#' @export varmatCUBE
#' @param ordinal Vector of ordinal responses 
#' @param m Number of ordinal categories
#' @param param Vector of parameters for the specified CUBE model
#' @param Y Matrix of selected covariates to explain the uncertainty component (default: no covariate is included 
#' in the model)
#' @param W Matrix of selected covariates to explain the feeling component (default: no covariate is included 
#' in the model)
#' @param Z Matrix of selected covariates to explain the overdispersion component (default: no covariate is included 
#' in the model)
#' @param expinform Logical: if TRUE  and no covariate is included in the model, the function returns
#'  the expected variance-covariance matrix (default is FALSE: the function returns the observed 
#'  variance-covariance matrix)
#' @details The function checks if the variance-covariance matrix is positive-definite: if not, 
#' it returns a warning message and produces a matrix with NA entries.  No missing value should be present neither
#'   for \code{ordinal} nor for covariate matrices: thus, deletion or imputation procedures should be preliminarily run.
#' @seealso  \code{\link{vcov}}, \code{\link{cormat}}
#' @keywords htest
#' @references Iannario, M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data, 
#' \emph{Communications in Statistics - Theory and Methods}, \bold{43}, 771--786 \cr
#' Piccolo D. (2015). Inferential issues for CUBE models with covariates,
#'  \emph{Communications in Statistics. Theory and Methods}, \bold{44}(23), 771--786. \cr
#' @examples
#' m<-7; n<-500
#' pai<-0.83; csi<-0.19; phi<-0.045
#' ordinal<-simcube(n,m,pai,csi,phi)
#' param<-c(pai,csi,phi)
#' varmat<-varmatCUBE(ordinal,m,param)
#' ##########################
#' ### Including covariates
#' data(relgoods)
#' m<-10
#' naord<-which(is.na(relgoods$Tv))
#' nacov<-which(is.na(relgoods$BirthYear))
#' na<-union(naord,nacov)
#' age<-2014-relgoods$BirthYear[-na]
#' lage<-log(age)-mean(log(age))
#' Y<-W<-Z<-lage
#' ordinal<-relgoods$Tv[-na]
#' estbet<-c(0.18,1.03); estgama<-c(-0.6,-0.3); estalpha<-c(-2.3,0.92)
#' param<-c(estbet,estgama,estalpha)
#' varmat<-varmatCUBE(ordinal,m,param,Y=Y,W=W,Z=Z,expinform=TRUE)





varmatCUBE<-function(ordinal,m,param,Y=0,W=0,Z=0,expinform=FALSE){

  if (is.factor(ordinal)){
    ordinal<-unclass(ordinal)
  }
  
  ry<-NROW(Y); rw<-NROW(W); rz<-NROW(Z);
  
  
  if(ry==1 & rw==1 & rz==1) {
    pai<-param[1]; csi<-param[2]; phi<-param[3];
    if(expinform==FALSE){
      freq<-tabulate(ordinal,nbins=m)
      varmat<-varcovcubeobs(m,pai,csi,phi,freq)
    } else{
      n<-length(ordinal)
      varmat<-varcovcubeexp(m,pai,csi,phi,n)
    }
  } else {
    if(ry>1 & rz>1 & rw >1){
      ncy<-NCOL(Y)
      ncw<-NCOL(W)
      Y<-as.matrix(Y);W<-as.matrix(W); Z<-as.matrix(Z);
      
      if (ncol(Y)==1){
        Y<-as.numeric(Y)
      }
      if (ncol(W)==1){
        W<-as.numeric(W)
      }
      if (ncol(Z)==1){
        Z<-as.numeric(Z)
      }
      
      estbet<-param[1:(ncy+1)]; estgama<-param[(ncy+2):(ncy+ncw+2)]; estalpha<-param[(ncy+ncw+3):length(param)];
      varmat<-varcovcubecov(m,ordinal,Y,W,Z,estbet,estgama,estalpha) 
    } else {
      cat("CUBE models not available for this variables specification")
    }
  }
  return(varmat)
  
}