#' @title Probability distribution of a CUSH model
#' @aliases probcush
#' @description Compute the probability distribution of a CUSH model without covariates, that is a mixture of a 
#' degenerate random variable with mass at the shelter category and the Uniform distribution.
#' @keywords distribution
#' @export probcush 
#' @usage probcush(m,delta,shelter)
#' @param m Number of ordinal categories
#' @param delta Shelter parameter
#' @param shelter Category corresponding to the shelter choice
#' @return The vector of the probability distribution of a CUSH model without covariates.
#' @references 
#' Capecchi S. and Piccolo D. (2017). Dealing with heterogeneity in ordinal responses,
#'  \emph{Quality and Quantity}, \bold{51}(5), 2375--2393 \cr
#' Capecchi S. and Iannario M. (2016). Gini heterogeneity index for detecting uncertainty in ordinal data surveys,
#'  \emph{Metron}, \bold{74}(2), 223--232
#' @examples
#' m<-10
#' shelter<-1
#' delta<-0.4
#' pr<-probcush(m,delta,shelter)
#' plot(1:m,pr,type="h",xlab="Number of categories")
#' points(1:m,pr,pch=19)



probcush <-
function(m,delta,shelter){
  delta*(ifelse(seq(1,m)==shelter,1,0) - 1/m) + 1/m
}
