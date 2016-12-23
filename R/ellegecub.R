#' @title Log-likelihood function for gecub distribution
#' @description Log-likelihood function for gecub distribution
#' @aliases ellegecub
#' @usage ellegecub(ordinal,Y,W,X,bet,gama,omega,shelter)
#' @param ordinal Vector of ordinal responses
#' @param Y Matrix of selected covariates for explaining the uncertainty component, not including intercept
#' @param W Matrix of selected covariates for explaining the feeling component, not including intercept
#' @param X Matrix of selected covariates for explaining the shelter effect, not including intercept
#' @param bet Matrix of selected covariates for explaining the uncertainty component, not including intercept
#' @param gama Matrix of selected covariates for explaining the feeling component, not including intercept
#' @param omega Matrix of selected covariates for explaining the shelter effect, not including intercept
#' @param shelter Category corresponding to the shelter choice
#' @keywords internal 


ellegecub<-function(ordinal,Y,W,X,bet,gama,omega,shelter){
  
  
  
  probn<-probgecub(factor(ordinal,ordered=TRUE),Y,W,X,bet,gama,omega,shelter);
  return(sum(log(probn)));
}
