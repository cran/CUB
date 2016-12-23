#' @title Probability distribution of an IHG model with covariates
#' @aliases probihgcovn
#' @description Given a vector of \eqn{n} ratings over \eqn{m} categories, it returns a vector
#'  of length \eqn{n} whose $i$-th element is the  probability of observing the $i$-th rating for the 
#'  corresponding IHG model with parameter \eqn{\theta_i}, obtained via logistic link with covariates
#'   and coefficients.
#' @keywords distribution
#' @export probihgcovn
#' @usage probihgcovn(m,ordinal,U,nu)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses (of factory type)
#' @param U Matrix of selected covariates for explaining the preference parameter
#' @param nu Vector of coefficients for covariates, whose length equals NCOL(U)+1 to include
#'  an intercept term in the model (first entry)
#' @details The matrix \eqn{U} is expanded with a vector with entries equal to 1 in the first column to include
#'  an intercept term in the model.
#' @seealso \code{\link{probihg}}
#' @examples
#' n<-100
#' m<-7
#' theta<-0.30
#' ordinal<-simihg(n,m,theta)
#' U<-sample(c(0,1),n,replace=TRUE)
#' nu<-c(0.12,-0.5)
#' pr<-probihgcovn(m,ordinal,U,nu)



probihgcovn <-
function(m,ordinal,U,nu){
  
  if (!is.factor(ordinal)){
    stop("Response must be an ordered factor")
  }
  
  ordinal<-unclass(ordinal)
  
  n<-length(ordinal)
  vett<-rep(NA,n)
  U<-as.matrix(U)
  if (ncol(U)==1){
    U<-as.numeric(U)
  }
 
  thetavett<-logis(U,nu)
  for (i in 1:n){
    prob<-probihg(m,thetavett[i])
    vett[i]<-prob[ordinal[i]]
  }
  return(vett)
}
