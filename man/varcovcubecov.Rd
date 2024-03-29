% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/varcovcubecov.R
\name{varcovcubecov}
\alias{varcovcubecov}
\title{Variance-covariance matrix of a CUBE model with covariates}
\usage{
varcovcubecov(m, ordinal, Y, W, Z, estbet, estgama, estalpha)
}
\arguments{
\item{m}{Number of ordinal categories}

\item{ordinal}{Vector of ordinal responses}

\item{Y}{Matrix of covariates for explaining the uncertainty component}

\item{W}{Matrix of covariates for explaining the feeling component}

\item{Z}{Matrix of covariates for explaining the overdispersion component}

\item{estbet}{Vector of the estimated parameters for the uncertainty component, with length equal to
NCOL(Y)+1 to account for an intercept term (first entry)}

\item{estgama}{Vector of the estimated parameters for the  feeling component, with length equal to
NCOL(W)+1 to account for an intercept term (first entry)}

\item{estalpha}{Vector of the estimated parameters for the overdispersion component, with length
equal to NCOL(Z)+1 to account for an intercept term (first entry)}
}
\description{
Compute the variance-covariance matrix of parameter estimates of a CUBE model with covariates
for all the three parameters.
}
\details{
The function checks if the variance-covariance matrix is positive-definite: if not,
it returns a warning message and produces a matrix with NA entries.
}
\references{
Piccolo, D. (2014), Inferential issues on CUBE models with covariates,
\emph{Communications in Statistics - Theory and Methods}, \bold{44}, DOI: 10.1080/03610926.2013.821487
}
\keyword{internal}
