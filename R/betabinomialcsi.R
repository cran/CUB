#' @title Beta-Binomial probabilities of ordinal responses, given feeling parameter for each observation
#' @description Compute the Beta-Binomial probabilities of the given ordinal responses, with feeling 
#' parameter specified for each observation, 
#' and with the same overdispersion parameter for all the responses.
#' @aliases betabinomialcsi
#' @usage betabinomialcsi(m,ordinal,csivett,phi)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses (factor type)
#' @param csivett  Vector of feeling parameters of the Beta-Binomial distribution for the given ordinal 
#' responses
#' @param phi Overdispersion parameter of the Beta-Binomial distribution 
#' @export betabinomialcsi
#' @return A vector of the same length as ordinal, containing the Beta-Binomial probability of each observation 
#' for the corresponding feeling parameter and for the specified overdispersion parameter
#' @seealso   \code{\link{betar}}, \code{\link{betabinomial}}
#' @references Iannario, M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data,
#'  \emph{Communications in Statistics - Theory and Methods}, \bold{43}, 771--786 \cr
#'  Piccolo D. (2015). Inferential issues for CUBE models with covariates. 
#' \emph{Communications in Statistics - Theory and Methods}, \bold{44}(23), 771--786.
#' @keywords distribution
#' @examples 
#' data(relgoods)
#' attach(relgoods)
#' m<-10
#' naTv<-which(is.na(Tv))
#' ordinal<-na.omit(Tv)
#' age<-2014-BirthYear[-naTv]
#' lage<-log(age)-mean(log(age))
#' W<-lage
#' gama<-c(-0.61,-0.31)
#' phi<-0.16 
#' csivett<-logis(W,gama)
#' pr<-betabinomialcsi(m,ordinal,csivett,phi)


betabinomialcsi <-function(m,ordinal,csivett,phi){
  
  if (!is.factor(ordinal)){
    stop("Response must be an ordered factor")
  }
  ordinal<-unclass(ordinal)
  
  n<-length(ordinal)
  betabin<-rep(NA,m)
  for(i in 1:n){
    bebeta<-betar(m,csivett[i],phi)
    betabin[i]<-bebeta[ordinal[i]]   
  }
  return(betabin)
}


