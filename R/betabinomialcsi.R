#' @title Beta-Binomial probabilities of ordinal responses, given feeling parameter for each observation
#' @description Compute the Beta-Binomial probabilities of given ordinal responses, with feeling 
#' parameter specified for each observation, 
#' and with the same overdispersion parameter for all the responses.
#' @aliases betabinomialcsi
#' @usage betabinomialcsi(m,ordinal,csivett,phi)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses. Missing values are not allowed: they should be preliminarily deleted 
#' or imputed
#' @param csivett  Vector of feeling parameters of the Beta-Binomial distribution for given ordinal 
#' responses
#' @param phi Overdispersion parameter of the Beta-Binomial distribution 
#' @export betabinomialcsi
#' @return A vector of the same length as ordinal: each entry is the Beta-Binomial probability for the given observation 
#' for the corresponding feeling and overdispersion parameters.
#' @seealso   \code{\link{betar}}, \code{\link{betabinomial}}
#' @references Iannario, M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data,
#'  \emph{Communications in Statistics - Theory and Methods}, \bold{43}, 771--786 \cr
#'  Piccolo D. (2015). Inferential issues for CUBE models with covariates. 
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
#' gama<-c(-0.61,-0.31)
#' phi<-0.16 
#' csivett<-logis(lage,gama)
#' pr<-betabinomialcsi(m,ordinal,csivett,phi)
#' plot(density(pr))


betabinomialcsi <-function(m,ordinal,csivett,phi){

  if (is.factor(ordinal)){
    ordinal<-unclass(ordinal)
  }
  n<-length(ordinal)
  betabin<-rep(NA,m)
  for(i in 1:n){
    bebeta<-betar(m,csivett[i],phi)
    betabin[i]<-bebeta[ordinal[i]]   
  }
  return(betabin)
}


