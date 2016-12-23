#' @title Pearson \eqn{X^2} statistic
#' @description Compute the \eqn{X^2} statistic of Pearson for CUB models with one or two discrete
#'  covariates for the feeling component.
#' @aliases chi2cub
#' @usage chi2cub(m,ordinal,W,pai,gama)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses (factor type)
#' @param W Matrix of covariates for the feeling component
#' @param pai Uncertainty parameter 
#' @param gama Vector of parameters for the feeling component, with length equal to NCOL(W)+1 
#' to account for an intercept term (first entry of \code{gama}) 
#' @export chi2cub
#' @return A list with the following components: 
#' \item{df}{Degrees of freedom}
#' \item{chi2}{Value of the Pearson fitting measure}
#' \item{dev}{Deviance indicator}
#' @keywords htest
#' @references 
#' Tutz, G. (2012). \emph{Regression for Categorical Data}, Cambridge University Press, Cambridge 
#' @examples 
#' data(univer)
#' attach(univer)
#' m<-7
#' pai<-0.3
#' gama<-c(0.1,0.7)
#' pearson<-chi2cub(m,informat,W=gender,pai,gama)
#' degfree<-pearson$df
#' statvalue<-pearson$chi2
#' deviance<-pearson$dev


chi2cub <-
function(m,ordinal,W,pai,gama){
  
  if (!is.factor(ordinal)){
    stop("Response must be an ordered factor")
  }
  
  W <- as.matrix(W)
  nw<-NCOL(W)
  if(nw==1){
    chi2cub1cov(m,ordinal,W,pai,gama)
  } else if(nw==2) {
    chi2cub2cov(m,ordinal,W[,1],W[,2],pai,gama)
  } else{
    cat("Works only for at most two covariates")
  }
  
}
