#' @title Pearson \eqn{X^2} statistic  for CUB models with one discrete covariate for feeling
#' @description Compute the \eqn{X^2} statistic of Pearson for the goodness of fit of a CUB model for ordinal responses, where the feeling parameter 
#' is explained via a logistic transform of the only discrete covariate. It groups ratings in 
#' classes according to the values of the covariate. 
#' @aliases chi2cub1cov
#' @usage chi2cub1cov(m, ordinal, covar, pai, gama)
#' @param m Integer: number of ordinal categories
#' @param ordinal Vector of ordinal responses 
#' @param covar Vector of the selected covariate for explaining the feeling component
#' @param pai Uncertainty parameter 
#' @param gama gama}{Vector of parameters for the feeling component, with length equal to 2
#'  to account for an intercept term (first entry)
#' @return It returns the following results in a list: 
#' \item{df}{Number of degrees of freedom}
#' \item{chi2}{Value of the Pearson fitting measure}
#' \item{dev}{Deviance indicator}
#' @keywords internal
#' @import stats
#' @references Tutz, G. (2011). \emph{Regression for categorical data}, Cambridge Series in Statistical 
#' and Probabilistic Mathematics



chi2cub1cov <-function(m,ordinal,covar,pai,gama){
  
  
  covar<-as.matrix(covar)
  
  n<-length(ordinal)
  elle<-as.numeric(sort(unique(covar)))
  
  kappa<-length(elle)
  
  matfrel<-matrix(NA,nrow=kappa,ncol=m)
  matprob<-matrix(NA,nrow=kappa,ncol=m)
  
  chi2<-0
  dev<-0
  
  j<-1
  while(j<=kappa){
    quali<-which(covar==elle[j])
    Wquali<-covar[quali]
    qualiord<-ordinal[quali]
    nk<- length(qualiord)
    matfrel[j,]=tabulate(qualiord,nbins=m)/nk
    nonzero<-which(matfrel[j,]!=0)
    paij<-pai
    csij<-1/(1+ exp(-gama[1]-gama[2]*elle[j]))
    matprob[j,]<-t(probcub00(m,paij,csij))
    chi2<-chi2+nk*sum(((matfrel[j,]-matprob[j,])^2)/matprob[j,])
    dev<- dev + 2*nk*sum(matfrel[j,nonzero]*log(matfrel[j,nonzero]/matprob[j,nonzero]))
    j<-j+1
  }
    
  df<- kappa*(m-1)-(length(gama)+1)
  cat("Degrees of freedom         ==>  df  =",df, "\n")
  cat("Pearson Fitting measure    ==>  X^2 =",chi2,"(p-val.=",1-pchisq(chi2,df),")","\n")
  cat("Deviance                   ==>  Dev =",dev,"(p-val.=",1-pchisq(dev,df),")","\n")
  results<-list('chi2'=chi2,'df'=df,'dev'=dev)
  
}
