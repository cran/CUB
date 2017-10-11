#' @title Logarithmic score
#' @description Compute the logarithmic score of a CUB model with covariates both for the uncertainty 
#' and the feeling parameters.
#' @aliases logscore
#' @export logscore
#' @usage logscore(m,ordinal,Y,W,bet,gama)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses 
#' @param Y Matrix of covariates for explaining the uncertainty component
#' @param W Matrix of covariates for explaining the feeling component
#' @param bet Vector of parameters for the uncertainty component, with length NCOL(Y)+1 
#' to account for an intercept term (first entry of \code{bet})
#' @param gama Vector of parameters for the feeling component, with length NCOL(W)+1 
#' to account for an intercept term (first entry of \code{gama})
#' @details No missing value should be present neither
#'   for \code{ordinal} nor for covariate matrices: thus, deletion or imputation procedures should be
#'    preliminarily run.
#' @references 
#' Tutz, G. (2012). \emph{Regression for Categorical Data}, Cambridge University Press, Cambridge
#' @keywords htest
#' @examples
#' data(relgoods)
#' m<-10
#' naord<-which(is.na(relgoods$Walking))
#' nacovpai<-which(is.na(relgoods$Gender))
#' nacovcsi<-which(is.na(relgoods$Smoking))
#' na<-union(naord,union(nacovpai,nacovcsi))
#' ordinal<-relgoods$Walking[-na]
#' Y<-relgoods$Gender[-na]
#' W<-relgoods$Smoking[-na]
#' bet<-c(-0.45,-0.48)
#' gama<-c(-0.55,-0.43)
#' logscore(m,ordinal,Y=Y,W=W,bet,gama)


logscore <-
function(m,ordinal,Y,W,bet,gama){
  

  
  pr<-as.numeric(probcubpq(m,ordinal,Y,W,bet,gama))
  return(-2*sum(log(pr)))
}
