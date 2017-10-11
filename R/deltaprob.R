#' @title Mean difference of a discrete random variable
#' @description Compute the Gini mean difference of a discrete distribution
#' @usage deltaprob(prob)
#' @aliases deltaprob
#' @param prob Vector of the probability distribution 
#' @keywords univar
#' @return Numeric value of the Gini mean difference of the input probability distribution,
#'  computed according to the de Finetti-Paciello formulation.
#' @export deltaprob
#' @examples
#' prob<-c(0.04,0.04,0.05,0.10,0.21,0.32,0.24)
#' deltaprob(prob)


deltaprob <-
function(prob){
  frip<-cumsum(prob)
  m<-length(prob)
  frip1<-frip[1:(m-1)]
  2*sum(frip1*(1-frip1))
}
