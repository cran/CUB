#' @title Auxiliary matrix
#' @description Returns an auxiliary matrix needed for computing the variance-covariance matrix of a CUBE model with covariates.
#' @aliases auxmat
#' @usage auxmat(m, vettcsi, vettphi, a, b, c, d, e)
#' @param m Number of ordinal categories
#' @param vettcsi Vector of the feeling parameters of the Beta-Binomial distribution, with length equal to the number of observations
#' @param vettphi Vector of the overdispersion parameters of the Beta-Binomial distribution, with length equal to the number of observations
#' @param a Real number
#' @param b Real number
#' @param c Real number
#' @param d Real number
#' @param e Real number
#' @keywords internal
#' @references
#'  Iannario, M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data,
#'   \emph{ Communications in Statistics- Theory and Methods}, \bold{43}, 771--786 \cr
#'  Piccolo, D. (2014). Inferential issues on CUBE models with covariates, 
#'  \emph{Communications in Statistics. Theory and Methods}, \bold{44}, DOI: 10.1080/03610926.2013.821487 


auxmat <-
function(m,vettcsi,vettphi,a,b,c,d,e){
  elemat<-matrix(NA,nrow=m,ncol=length(vettcsi))
  for(k in 1:m){
    elemat[k,]<-e*((k-1)^d)/((a+b*vettcsi+vettphi*(k-1))^c)
  }
  return(elemat)
}


