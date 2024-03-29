% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitted.R
\name{fitted.GEM}
\alias{fitted.GEM}
\title{S3 method "fitted" for class "GEM"}
\usage{
\method{fitted}{GEM}(object, ...)
}
\arguments{
\item{object}{An object of class \code{\link{GEM}}}

\item{...}{Other arguments}
}
\description{
S3 method fitted for objects of class \code{\link{GEM}}.
}
\details{
Returns the fitted probability distribution for GEM models with no covariates. If only one dichotomous
covariate is included in the model to explain some components, it returns the fitted probability distribution for each profile.
}
\examples{
fitcub<-GEM(Formula(global~0|freqserv|0),family="cub",data=univer)
fitted(fitcub,digits=4)
}
\seealso{
\code{GEM}
}
\keyword{package}
