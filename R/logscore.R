#' @title Logarithmic score
#' @description Compute the logarithmic score of a CUB model with covariates both for the uncertainty 
#' and the feeling parameters.
#' @aliases logscore
#' @export logscore
#' @usage logscore(m,ordinal,Y,W,bet,gama)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses (factor type)
#' @param Y Matrix of covariates for explaining the uncertainty component
#' @param W Matrix of covariates for explaining the feeling component
#' @param bet Vector of parameters for the uncertainty component, with length NCOL(Y)+1 
#' to account for an intercept term (first entry of \code{bet})
#' @param gama Vector of parameters for the feeling component, with length NCOL(W)+1 
#' to account for an intercept term (first entry of \code{gama})
#' @references 
#' Tutz, G. (2012). \emph{Regression for Categorical Data}, Cambridge University Press, Cambridge
#' @keywords htest
#' @examples
#' data(relgoods)
#' attach(relgoods)
#' m<-10
#' naord<-which(is.na(Walking))
#' nacovpai<-which(is.na(Gender))
#' nacovcsi<-which(is.na(Smoking))
#' na<-union(naord,union(nacovpai,nacovcsi))
#' ordinal<-Walking[-na]
#' Y<-Gender[-na]
#' W<-Smoking[-na]
#' bet<-c(-0.45,-0.48)
#' gama<-c(-0.55,-0.43)
#' logscore(m,ordinal,Y=Y,W=W,bet,gama)


logscore <-
function(m,ordinal,Y,W,bet,gama){
  
  if (!is.factor(ordinal)){
    stop("Response must be an ordered factor")
  }
  
  pr<-as.numeric(probcubpq(m,ordinal,Y,W,bet,gama))
  return(-2*sum(log(pr)))
}
