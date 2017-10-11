#' @title Normalized Gini heterogeneity index
#' @description Compute the normalized Gini heterogeneity 
#' index for a given discrete probability distribution.
#' @usage gini(prob)
#' @aliases gini
#' @param prob Vector of  probability distribution or relative frequencies
#' @keywords univar
#' @export gini
#' @seealso \code{\link{laakso}}
#' @examples 
#' prob<-c(0.04,0.04,0.05,0.10,0.21,0.32,0.24)
#' gini(prob)


gini <-
  function(prob){
    m<-length(prob)
    (1-sum(prob^2))*m/(m-1)
  }
