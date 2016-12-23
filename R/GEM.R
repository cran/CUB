#' @title Main function for GEM models
#' @description Main function to estimate and validate Generalized Mixture models with uncertainty. 
#' @aliases GEM
#' @usage GEM(Formula,family=c("cub","cube","ihg","cush"),data,...)
#' @param Formula Object of class Formula. Response variable is the vector of ordinal responses - see Details.
#' @param family Character string indicating which class of GEM models to fit.
#' @param data an optional data frame (or object coercible by \code{as.data.frame} to a data frame)
#'  containing the variables in the model. If not found in data, the variables are taken from \code{environment(Formula)}.
#' @param ... Additional arguments to be passed for the specification of the model. See details.
#' @export GEM
#' @return An object of the class "GEM" is a list containing the following elements: 
#' \item{estimates}{Maximum likelihood estimates of parameters}
#' \item{loglik}{Log-likelihood function at the final estimates}
#' \item{varmat}{Variance-covariance matrix of final estimates}
#' \item{niter}{Number of executed iterations}
#' \item{BIC}{BIC index for the estimated model}
#' \item{ordinal}{Vector of ordinal responses on which the model has been fitted}
#' \item{time}{Processor time for execution}
#' \item{ellipsis}{Retrieve the arguments passed to the call and extra arguments generated via the call}
#' \item{family}{Character string indicating the sub-class of the fitted model}
#' \item{formula}{Returns the Formula of the call for the fitted model}
#' \item{call}{Returns the executed call}
#' @import methods Formula
#' @details It is the main function for GEM models, calling for the corresponding function for
#'  the specified subclass. The number of category is obtained as the number of levels of the ordinal observations. 
#'  If  \code{family="cub"}, then a CUB mixture model is fitted to the data to explain uncertainty, 
#'  feeling and possible shelter effect by further passing the extra argument \code{shelter} for the corresponding category.
#'  Subjects' covariates can be included by specifying covariates matrices in the 
#'  Formula as \code{ordinal~Y|W|X},  to explain uncertainty (Y), feeling (W) or shelter (X).  Notice that
#'  covariates for shelter effect can be included only if specified for both feeling and uncertainty. \cr
#'  If \code{family="cube"}, then a CUBE mixture model is fitted to the data to explain uncertainty, feeling and overdispersion.
#'  Subjects' covariates can be also included to explain the feeling component or all the three components by
#'  specifying covariates matrices in the Formula as \code{ordinal~Y|W|Z} to explain uncertainty (Y), feeling (W) or 
#'  overdispersion (Z). An extra logical argument \code{expinform} indicates whether or not to use the expected or the 
#'  observed information matrix (default is FALSE). \cr
#'  If \code{family="ihg"}, then an IHG model is fitted to the data. IHG models are nested into CUBE models
#'  (see the references below). The parameter \eqn{\theta} gives the probability of observing a rating corresponding 
#'  to the first category and is therefore a direct measure of preference, attraction, pleasantness toward the 
#'  investigated item. This is the reason why \eqn{\theta} is customarily referred to as the preference parameter of the 
#'  IHG model. Covariates for the preference parameter \eqn{\theta} have to be specified in matrix form in the Formula as \code{ordinal~U}. 
#'  If \code{family="cush"}, then a CUSH model is fitted to the data. The category corresponding to the inflation should be
#'  passed via argument shelter=shelter. Covariates for the shelter parameter \eqn{\delta}
#'  are specified in matrix form Formula as \code{ordinal~X}. \cr
#'  Even if no covariate is included in the model for a given component, the corresponding model matrix needs always
#'  to be specified: in this case, it should be set to 0 (see examples below). Extra arguments include the maximum 
#'  number of iterations (\code{maxiter}, default: \code{maxiter}=500) for the optimization algorithm and 
#'  the required error tolerance (\code{toler}, default: \code{toler}=1e-6). \cr
#'  Standard methods: \code{logLik()}, \code{BIC()}, \code{vcov()}, \code{fitted()}, \code{coef()}, \code{print()}, \code{summary()}
#'  are implemented.\cr
#'  The optimization procedure is run via \code{optim()} when required. If the estimated variance-covariance matrix is not positive definite, the function returns a 
#'  warning message and produces a matrix with NA entries.
#' @references 
#' Iannario M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data, 
#' \emph{Communications in Statistics - Theory and Methods}, \bold{43}, 771--786 \cr
#' Piccolo D. (2015). Inferential issues for CUBE models with covariates,
#'  \emph{Communications in Statistics. Theory and Methods}, \bold{44}(23), 771--786. \cr
#' Iannario M. (2015). Detecting latent components in ordinal data with overdispersion by means
#'  of a mixture distribution, \emph{Quality & Quantity}, \bold{49}, 977--987 \cr
#'  Iannario M. (2016). Testing the overdispersion parameter in CUBE models. 
#'  \emph{Communications in Statistics: Simulation and Computation}, \bold{45}(5), 1621--1635.\cr
#'  Iannario M. and Piccolo D. (2016a). A comprehensive framework for regression models of ordinal data.
#'  \emph{Metron}, \bold{74}(2), 233--252.\cr
#'  Iannario M. and Piccolo D. (2016b). A generalized framework for modelling ordinal data. 
#'  \emph{Statistical Methods and Applications}, \bold{25}, 163--189.\cr 
#'  D'Elia A. (2003). Modelling ranks using the inverse hypergeometric distribution, 
#' \emph{Statistical Modelling: an International Journal}, \bold{3}, 65--78 \cr
#' Capecchi S. and Piccolo D. (2016). Dealing with heterogeneity in ordinal responses,
#'  \emph{Quality and Quantity}, DOI: 10.1007/s11135-016-0393-3 \cr
#' @seealso \code{\link{logLik}}, \code{\link{coef}}, \code{\link{BIC}},  \code{\link{makeplot}},
#'  \code{\link{summary}}, \code{\link{vcov}}, \code{\link{fitted}}, \code{\link{cormat}}
#' @keywords models
#' @examples 
#' \donttest{
#' ## CUB models with no covariates
#' data(relgoods)
#' attach(relgoods)
#' F1<-Formula(Walking~0|0|0)
#' model<-GEM(F1,family="cub")     
#' coef(model,digits=5)     # Estimated parameter vector (pai,csi)
#' logLik(model)            # Log-likelihood function at ML estimates
#' vcov(model,digits=4)     # Estimated Variance-Covariance matrix
#' cormat(model)            # Parameter Correlation matrix
#' fitted(model)            # Fitted probability distribution
#' nniter<-model$niter  
#' ################
#' ## CUB model with shelter effect
#' data(univer)
#' attach(univer)
#' F2<-Formula(officeho~0|0|0)
#' model<-GEM(F2,family="cub",shelter=7)
#' BICcub<-BIC(model,digits=4)
#' ################
#' ## CUB model with covariate for all components - GeCUB
#' dichoage<-ifelse(age<=25,0,1)
#' F3<-Formula(officeho~gender|dichoage|freqserv)
#' modelgecub<-GEM(F3,family="cub",shelter=7)
#' BICgecub<-BIC(modelgecub)
#' ################
#' ## CUB model with covariate for uncertainty
#' data(relgoods)
#' attach(relgoods)
#' naord<-which(is.na(Parents))
#' nacov<-which(is.na(Smoking))
#' na<-union(naord,nacov)
#' ordinal<-Parents[-na]
#' cov<-Smoking[-na]
#' modelcovpai<-GEM(Formula(ordinal~cov|0|0),family="cub")
#' fitted(modelcovpai)
#' makeplot(modelcovpai)
################
#' ## CUB model with covariate for feeling
#' data(relgoods)
#' attach(relgoods)
#' naord<-which(is.na(RelFriends))
#' nacov<-which(is.na(WalkAlone))
#' na<-union(naord,nacov)
#' ordinal<-RelFriends[-na]
#' cov<-WalkAlone[-na]
#' modelcovcsi<-GEM(Formula(ordinal~0|cov|0),family="cub")
#' fitted(modelcovcsi)
#' makeplot(modelcovcsi)
#' # all methods 
#' ##################
#' ## CUB model with covariates for both uncertainty and feeling components
#' data(univer)
#' attach(univer)
#' lage<-log(age)-mean(log(age))
#' F6<-Formula(global~gender|lage|0)
#' model<-GEM(Formula=F6,family="cub",data=univer,maxiter=150,toler=1e-4) 
#' param<-coef(model)
#' bet<-param[1:2]      # ML estimates of coefficients for uncertainty covariate: gender
#' gama<-param[3:4]     # ML estimates of coefficients for feeling covariate: lage
#' ################################################ 
#' # CUBE models with no covariates
#' data(relgoods)
#' attach(relgoods)
#' naord<-which(is.na(MeetRelatives))   
#' nacov<-which(is.na(WalkAlone))
#' na<-union(naord,nacov)
#' ordinal<-MeetRelatives[-na]
#' cov<-WalkAlone[-na]
#' m<-length(levels(ordinal))
#' starting<-inibestcube(m,ordinal)
#' Fcube000<-Formula(MeetRelatives~0|0|0)
#' model<-GEM(Fcube000,family="cube",starting=starting,expinform=TRUE) 
#' coef(model,digits=4)       # Final ML estimates
#' logLik(model)     # Maximum value of the log-likelihood function
#' vcov(model)           
#' print(model)
#' BIC(model)
#' fitted(model)
#' makeplot(model)
#' summary(model)
#' ##########
#' ## CUBE with covariate WalkAlone only for the feeling component
#' Fcube010<-Formula(ordinal~0|cov|0)
#' modelcovcsi<-GEM(Fcube010,family="cube",maxiter=50)
#' summary(modelcovcsi)
#' ## CUBE with covariates for all components
#' Fcube111<-Formula(ordinal~cov|cov|cov)
#' modelcov<-GEM(Fcube111,family="cube",maxiter=50)
#' BIC(modelcovcsi)    #  modelcovcsi$BIC
#' BIC(modelcov)       #  modelcov$BIC
#' #######################################
#' # IHG models without covariates
#' data(univer)
#' attach(univer)
#' Fihg0<-Formula(willingn~0)
#' model<-GEM(Fihg0,family="ihg")
#' coef(model)                ## ML estimate of preference parameter theta
#' fitted(model)              ## fitted probabilities
#' makeplot(model)
#' #### with covariates
#' Fihgcov<-Formula(willingn~freqserv)
#' modelcov<-GEM(Formula=Fihgcov,family="ihg")
#' omega<-coef(modelcov)      ## ML estimates (intercept term: omega[1])  (modelcov$estimates)
#' maxlik<-logLik(modelcov)   ## also modelcov$loglik
#' makeplot(modelcov)
#' summary(modelcov)
#' #################################
#' # CUSH models without covariate
#' data(relgoods)
#' attach(relgoods)
#' Fcush0<-Formula(Dog~0)
#' model<-GEM(Fcush0,family="cush",shelter=1)
#' delta<-coef(model)      # ML estimates of delta
#' maxlik<-logLik(model)   # Log-likelihood at ML estimates
#' summary(model)
#' makeplot(model)
#' ###############################################
#' ### CUSH model with covariates
#' naord<-which(is.na(MusicInstr))
#' nacov<-which(is.na(Smoking))
#' na<-union(naord,nacov)
#' ordinal<-MusicInstr[-na]
#' cov<-Smoking[-na]
#' modelcov<-GEM(Formula(ordinal~cov),family="cush",shelter=1)
#' omega<-coef(modelcov)
#' maxlik<-logLik(modelcov)
#' varmat<-vcov(modelcov)
#' summary(modelcov,digits=3)
#' fitted(modelcov,digits=3)
#' }
#' 

GEM<-function(Formula,family=c("cub","cube","ihg","cush"),data,...){
  
  #... for maxiter, toler, shelter,Y, W, for cub etc
  # if (!is.factor(ordinal))  stop("Response must be an ordered factor")
  # lev <- levels(ordinal); m <- length(lev)
  # ordinal <- unclass(ordinal)
  
  call <- match.call()
  
  if (missing(data)) data <- environment(Formula)
  
  ellipsis.arg<-list(...)

  # if (!is.factor(ordinal))  stop("Response must be an ordered factor")

  mf<-model.frame(Formula,data)
  ordinal<-as.numeric(model.response(mf))
  
  lev <- levels(factor(ordinal,ordered=TRUE)); 
  m <- length(lev)  
  
  ellipsis.arg$m <-m
  
  #as.list(mapply(c, lista1, m, SIMPLIFY = TRUE))
  
  # # #cl<-match.call(expand.dots=TRUE)
  
  if (family == "cub"){
    if (is.null(ellipsis.arg$maxiter)){
      ellipsis.arg$maxiter <- 500
    }
    if (is.null(ellipsis.arg$toler)){
      ellipsis.arg$toler <- 1e-4
    }
    mod<-CUB(Formula,data,ellipsis.arg)
    modelname<-"CUB"
    
  }
  if (family =="cube"){
    if (is.null(ellipsis.arg$maxiter)){
      ellipsis.arg$maxiter <- 1000
    }
    if (is.null(ellipsis.arg$toler)){
      ellipsis.arg$toler <- 1e-6
    }
    maxiter<- ellipsis.arg$maxiter
    toler<-ellipsis.arg$toler
    mod<-CUBE(Formula,data,ellipsis.arg)
    modelname<-"CUBE"
    
  }
  if (family == "ihg"){
    mod<-IHG(Formula,data,...)
    modelname<-"IHG"
    ellipsis.arg$maxiter<-1
    ellipsis.arg$niter<-1
  }
  if (family == "cush"){
    mod<-CUSH(Formula,data,...)
    modelname<-"CUSH"
    ellipsis.arg$maxiter<-1
    ellipsis.arg$niter<-1
  }
  
  stime<-mod$estimates
  durata<-mod$time
  #maxiter<-mod$maxiter
  loglik<-as.numeric(mod$loglik)
  niter<-mod$niter
  varmat<-mod$varmat
  BIC<-as.numeric(mod$BIC)
  time<-mod$durata
  

  results<-list('estimates'=stime,'ordinal'=ordinal,'time'=durata,
                'loglik'=loglik,'niter'=niter,'varmat'=varmat,
                'BIC'=BIC,'ellipsis'=ellipsis.arg,'family'=modelname,
                'formula'=Formula,'call'=call)
  
  attr(results,"hidden")<-c("estimates","ordinal","loglik","varmat","BIC","ellipsis","family")
  
  class(results)<-"GEM"
  
  

  
  
  # print.GEM<-function(x,...){
  #   hid<-attr(x,"hidden")
  #   print(x[!names(x)%in%hid])
  # }
  # 
  

  if (family == "cub"){
    class(results)<-append("GEM","CUB")
    #class(results)<-c("CUB","GEM")
  }
  if (family =="cube"){
    class(results)<-append("GEM","CUBE")
  }
  if (family == "ihg"){
    class(results)<-append("GEM","IHG")
  }
  if (family == "cush"){
    class(results)<-append("GEM","CUSH")
  }
  
  invisible(results)
  #return(results)
}

