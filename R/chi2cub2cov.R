#' @title Pearson \eqn{X^2} statistic  for CUB models with one discrete covariate for feeling
#' @description Compute the \eqn{X^2} statistic of Pearson for the goodness of fit of a CUB model for ordinal responses, where the feeling parameter 
#' is explained via a logistic transform of the only discrete covariate. It groups ratings in 
#' classes according to the values of the covariate. 
#' @aliases chi2cub2cov
#' @usage chi2cub2cov(m, ordinal, covar1, covar2, pai, gama)
#' @param m Integer: number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param covar1 Vector of the first covariate values for explaining the feeling component
#' @param covar2 Vector of the second covariate values for explaining the feeling component
#' @param pai Uncertainty parameter 
#' @param gama Vector of parameters for the feeling component, with length equal to 2
#'  to account for an intercept term (first entry)
#' @return It returns the following results in a list: 
#' \item{df}{Number of degrees of freedom}
#' \item{chi2}{Value of the Pearson fitting measure}
#' \item{dev}{Deviance indicator}
#' @keywords internal
#' @import stats
#' @references
#' Tutz, G. (2011). \emph{Regression for categorical data}, Cambridge Series in Statistical 
#' and Probabilistic Mathematics


chi2cub2cov <-function(m,ordinal,covar1,covar2,pai,gama){
  
  
  n<-length(ordinal)
  W<-cbind(covar1,covar2)
  elle1<-as.numeric(sort(unique(covar1)))
  elle2<-as.numeric(sort(unique(covar2)))
  
  profiles<-expand.grid(elle1,elle2)
  kappa<-nrow(profiles)
  
  matfrel<-matrix(NA,nrow=kappa,ncol=m)
  matprob<-matrix(NA,nrow=kappa,ncol=m)
  
  chi2<-0
  dev<-0
  
  j<-1
  while(j<=kappa){
    quali<-which(covar1==profiles[j,1] & covar2==profiles[j,2])
    if (length(quali) == 0){
      j<-j+1; print(paste("Profile ",j-1,"is void")); cat("\n") ;
    }else {
      Wquali<-W[quali,]
      qualiord<-ordinal[quali]
      nk<- length(qualiord)
      matfrel[j,]<-tabulate(qualiord,nbins=m)/nk
      nonzero<-which(matfrel[j,]!=0)
      paij<-pai
      csij<-1/(1+ exp(-gama[1]-gama[2]*profiles[j,1]-gama[3]*profiles[j,2]))
      matprob[j,]<-t(probcub00(m,paij,csij))
      chi2<-chi2+nk*sum(((matfrel[j,]-matprob[j,])^2)/matprob[j,])
      dev<- dev + 2*nk*sum(matfrel[j,nonzero]*log(matfrel[j,nonzero]/matprob[j,nonzero]))
      j<-j+1
    }
    
  }
  
  df<-kappa*(m-1)-(length(gama)+1)
  cat("Degrees of freedom         ==>  df  =",df, "\n")
  cat("Pearson Fitting measure    ==>  X^2 =",chi2,"(p-val.=",1-pchisq(chi2,df),")","\n")
  cat("Deviance                   ==>  Dev =",dev,"(p-val.=",1-pchisq(dev,df),")","\n")
  results<-list('chi2'=chi2,'df'=df,'dev'=dev)
}
