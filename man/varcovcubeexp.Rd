% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/varcovcubeexp.R
\name{varcovcubeexp}
\alias{varcovcubeexp}
\title{Variance-covariance matrix for CUBE models based on the expected information matrix}
\usage{
varcovcubeexp(m, pai, csi, phi, n)
}
\arguments{
\item{m}{Number of ordinal categories}

\item{pai}{Uncertainty parameter}

\item{csi}{Feeling parameter}

\item{phi}{Overdispersion parameter}

\item{n}{Number of observations}
}
\description{
Compute the variance-covariance matrix of parameter estimates as the inverse of
the expected information matrix for a CUBE model without covariates.
}
\details{
The function checks if the variance-covariance matrix is positive-definite: if not,
it returns a warning message and produces a matrix with NA entries.
}
\references{
Iannario, M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data,
\emph{Communications in Statistics - Theory and Methods}, \bold{43}, 771--786
}
\seealso{
\code{\link{varcovcubeobs}}
}
\keyword{internal}
