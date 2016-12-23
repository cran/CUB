#' @title Variance-covariance matrix for CUB models
#' @aliases varmatCUB
#' @description Compute the variance-covariance matrix of parameter estimates for CUB models with or without
#' covariates for the feeling and the overdispersion parameter, and for extended CUB models with shelter effect.
#' @usage varmatCUB(ordinal,m,param,Y=0,W=0,X=0,shelter=0)
#' @export varmatCUB
#' @param ordinal Vector of ordinal responses (of factory type)
#' @param m Number of ordinal categories
#' @param param Vector of parameters for the specified CUB model
#' @param Y Matrix of selected covariates to explain the uncertainty component (default: no covariate is included 
#' in the model)
#' @param W Matrix of selected covariates to explain the feeling component (default: no covariate is included 
#' in the model)
#' @param X Matrix of selected covariates to explain the shelter effect (default: no covariate is included 
#' in the model)
#' @param shelter Category corresponding to the shelter choice (default: no shelter effect is included in the 
#' model)
#' @details The function checks if the variance-covariance matrix is positive-definite: if not, 
#' it returns a warning message and produces a matrix with NA entries.
#' @seealso  \code{\link{vcov}}, \code{\link{cormat}}
#' @keywords htest
#' @references Piccolo D. (2006). Observed Information Matrix for MUB Models, 
#' \emph{Quaderni di Statistica}, \bold{8}, 33--78 \cr
#' Iannario, M. (2012). Modelling shelter choices in ordinal data surveys. 
#' \emph{Statistical Modelling and Applications}, \bold{21}, 1--22 \cr
#' Iannario M. and Piccolo D. (2016b). A generalized framework for modelling ordinal data. 
#'  \emph{Statistical Methods and Applications}, \bold{25}, 163--189.\cr 
#' @examples
#' data(univer)
#' attach(univer)
#' m<-7
#' ### CUB model with no covariate
#' pai<-0.87; csi<-0.17 
#' param<-c(pai,csi)
#' varmat<-varmatCUB(global,m,param)
#' #######################
#' ### and with covariates for feeling
#' data(univer)
#' m<-7
#' pai<-0.86; gama<-c(-1.94,-0.17)
#' param<-c(pai,gama)
#' W<-gender      
#' varmat<-varmatCUB(willingn,m,param,W=W)
#' #######################
#' ### CUB model with uncertainty covariates
#' #' data(relgoods)
#' attach(relgoods)
#' m<-10
#' naord<-which(is.na(Physician))
#' nacov<-which(is.na(Gender))
#' na<-union(naord,nacov)
#' ordinal<-Physician[-na]
#' Y<-Gender[-na]
#' bet<-c(-0.81,0.93); csi<-0.20
#' varmat<-varmatCUB(ordinal,m,param=c(bet,csi),Y=Y)
#' #######################
#' ### and with covariates for both parameters
#' data(relgoods)
#' attach(relgoods)
#' m<-10
#' naord<-which(is.na(Physician))
#' nacov<-which(is.na(Gender))
#' na<-union(naord,nacov)
#' ordinal<-Physician[-na]
#' W<-Y<-Gender[-na]
#' gama<-c(-0.91,-0.7); bet<-c(-0.81,0.93)
#' varmat<-varmatCUB(ordinal,m,param=c(bet,gama),Y=Y,W=W)
#' #######################
#' ### Variance-covariance for a CUB model with shelter
#' m<-8; n<-300
#' pai1<-0.5; pai2<-0.3; csi<-0.4
#' shelter<-6
#' pr<-probcubshe1(m,pai1,pai2,csi,shelter)
#' ordinal<-factor(sample(1:m,n,prob=pr,replace=TRUE),ordered=TRUE)
#' param<-c(pai1,pai2,csi)
#' varmat<-varmatCUB(ordinal,m,param,shelter=shelter)







varmatCUB<-function(ordinal,m,param,Y=0,W=0,X=0,shelter=0){
  

  if (!is.factor(ordinal)){
    stop("Response must be an ordered factor")
  }
  
  ordinal<-unclass(ordinal)
  
  ry<-NROW(Y);   rw<-NROW(W); rx<-NROW(X); shelter<-as.numeric(shelter)
  
  if(shelter!=0){
    if (ry==1 & rw==1 & rx==1){
      pai1<-param[1]
      pai2<-param[2]
      csi<-param[3]
      n<-length(ordinal)
      varmat<-varcovcubshe(m,pai1,pai2,csi,shelter,n)
    }else if(ry!=1 & rw!=1 & rx!=1){
      ny<-NCOL(Y);nw<-NCOL(W);nx<-NCOL(X)
      Y<-as.matrix(Y);W<-as.matrix(W); X<-as.matrix(X);
      if (ncol(W)==1){
        W<-as.numeric(W)
      }
      if (ncol(Y)==1){
        Y<-as.numeric(Y)
      }
      if (ncol(X)==1){
        X<-as.numeric(X)
      }
      
      bet<-param[1:(ny+1)]; gama<-param[(ny+2):(ny+nw+2)]; omega<-param[(ny+nw+3):(ny+nw+nx+3)]
      varmat<-varcovgecub(ordinal,Y,W,X,bet,gama,omega,shelter)
      
    } else {
      varmat<-NULL
      cat("CUB model with shelter effect available only with covariates for all components")
    }
  }else{
    if(ry==1 & rw==1 & rx==1) {
      pai<-param[1]
      csi<-param[2]
      varmat<-varcovcub00(m,ordinal,pai,csi)
    }    else{
      if(ry!=1 & rw==1 & rx==1) {
        ncy<-NCOL(Y)
        Y<-as.matrix(Y);
        if (ncol(Y)==1){
          Y<-as.numeric(Y)
        }
        bet<-param[1:(ncy+1)]
        csi<-param[length(param)]
        varmat<-varcovcubp0(m,ordinal,Y,bet,csi)
      } else {
        if(ry==1 & rw!=1 & rx==1) {
          pai<-param[1]
          gama<-param[2:length(param)]
          W<-as.matrix(W);
          if (ncol(W)==1){
            W<-as.numeric(W)
          }
          
          varmat<-varcovcub0q(m,ordinal,W,pai,gama)
        } else{
          if(ry!=1 & rw!=1& rx==1) {
            Y<-as.matrix(Y);W<-as.matrix(W);ncy<-NCOL(Y)
            if (ncol(Y)==1){
              Y<-as.numeric(Y)
            }
            if (ncol(W)==1){
              W<-as.numeric(W)
            }
            bet<-param[1:(ncy+1)]
            gama<-param[(ncy+2):length(param)]
            varmat<-varcovcubpq(m,ordinal,Y,W,bet,gama)
          } else {
            cat("Wrong variables specification")
            varmat<-NULL
          }
        }                            
      }
    }
  }
  return(varmat)
}