#' @title Log-likelihood function for the IHG model with covariates
#' @aliases loglikihgcov
#' @description Compute the log-likelihood function for the IHG model
#'  with covariates to explain the preference parameter.
#' @usage loglikihgcov(m, ordinal, U, nu)
#' @keywords internal
#' @seealso loglikIHG
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param U Matrix of selected covariates for explaining the preference parameter
#' @param nu Vector of coefficients for covariates, whose length equals NCOL(U)+1 to include
#'  an intercept term in the model (first entry of nu)


loglikihgcov <-
function(m,ordinal,U,nu){
  
  if (is.factor(ordinal)){
    ordinal<-unclass(ordinal)
  }
  U<-as.matrix(U)
  sum(log(probihgcovn(m,ordinal,U,nu)))
}
