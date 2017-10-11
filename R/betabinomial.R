#' @title Beta-Binomial probabilities of ordinal responses, with feeling and overdispersion parameters
#' for each observation
#' @description Compute the Beta-Binomial probabilities of ordinal responses, given feeling and overdispersion
#' parameters for each observation.
#' @aliases betabinomial
#' @usage betabinomial(m,ordinal,csivett,phivett)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses. Missing values are not allowed: they should be preliminarily deleted 
#' or imputed
#' @param csivett  Vector of feeling parameters of the Beta-Binomial distribution for given ordinal responses
#' @param phivett Vector of overdispersion parameters of the Beta-Binomial distribution for given ordinal 
#' responses
#' @export betabinomial
#' @return A vector of the same length as ordinal, containing the Beta-Binomial probabilities of each observation,
#'  for the corresponding feeling and overdispersion parameters.
#' @details The Beta-Binomial distribution is the Binomial distribution in which the probability of success at
#'  each trial is random and follows the Beta distribution. It is frequently used in Bayesian 
#'  statistics, empirical Bayes methods and classical statistics as an overdispersed binomial distribution. 
#' @seealso  \code{\link{betar}}, \code{\link{betabinomialcsi}}
#' @references Iannario, M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data, 
#' \emph{Communications in Statistics - Theory and Methods}, \bold{43}, 771--786 \cr
#' Piccolo D. (2015). Inferential issues for CUBE models with covariates. 
#' \emph{Communications in Statistics - Theory and Methods}, \bold{44}(23), 771--786.
#' @keywords distribution
#' @examples 
#' data(relgoods)
#' m<-10
#' ordinal<-relgoods$Tv
#' age<-2014-relgoods$BirthYear
#' no_na<-na.omit(cbind(ordinal,age))
#' ordinal<-no_na[,1]; age<-no_na[,2]
#' lage<-log(age)-mean(log(age))
#' gama<-c(-0.6, -0.3)
#' csivett<-logis(lage,gama)
#' alpha<-c(-2.3,0.92); 
#' ZZ<-cbind(1,lage)
#' phivett<-exp(ZZ%*%alpha)
#' pr<-betabinomial(m,ordinal,csivett,phivett)
#' plot(density(pr))

betabinomial <-
function(m,ordinal,csivett,phivett){
  
  if (is.factor(ordinal)){
   ordinal<-unclass(ordinal)
  }
  
  n<-length(ordinal)
  betabin<-rep(NA,n)
  for(i in 1:n){
    bebeta<-betar(m,csivett[i],phivett[i])
    betabin[i]<-bebeta[ordinal[i]]   
  }
  return(betabin)
}

