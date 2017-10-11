#' @title Normalized Laakso and Taagepera heterogeneity index
#' @description Compute the normalized Laakso and Taagepera heterogeneity index for a given 
#' discrete probability distribution.
#' @aliases laakso
#' @export laakso
#' @usage laakso(prob)
#' @param prob Vector of a probability or relative frequency distribution
#' @seealso \code{\link{gini}}
#' @keywords univar
#' @references 
#' Laakso, M. and Taagepera, R. (1989). Effective number of parties: a measure with application to West Europe, 
#' \emph{Comparative Political Studies}, \bold{12}, 3--27.
#' @examples
#' prob<-c(0.04,0.04,0.05,0.10,0.21,0.32,0.24)
#' laakso(prob)



laakso <-
function(prob){
  m<-length(prob)
  1/(m/gini(prob)-m+1)}
