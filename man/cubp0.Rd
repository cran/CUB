% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cubp0.R
\name{cubp0}
\alias{cubp0}
\title{Main function for CUB models with covariates for the uncertainty component}
\usage{
cubp0(m, ordinal, Y, maxiter, toler)
}
\arguments{
\item{m}{Number of ordinal categories}

\item{ordinal}{Vector of ordinal responses}

\item{Y}{Matrix of selected covariates for explaining the uncertainty component}

\item{maxiter}{Maximum number of iterations allowed for running the optimization algorithm}

\item{toler}{Fixed error tolerance for final estimates}
}
\value{
An object of the class "CUB"
}
\description{
Estimate and validate a CUB model for given ordinal responses, with covariates for explaining
the feeling component via a logistic transform.
}
\references{
Iannario M. and Piccolo D. (2010), A new statistical model for the analysis of customer satisfaction,
\emph{Quality Technology and Quantity management}, \bold{7}(2) 149--168 \cr
Iannario M. and Piccolo D. (2012). CUB models: Statistical methods and empirical evidence, in:
Kenett R. S. and Salini S. (eds.), \emph{Modern Analysis of Customer Surveys: with applications using R},
J. Wiley and Sons, Chichester, 231--258
}
\keyword{internal}
