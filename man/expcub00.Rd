% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expcub00.R
\name{expcub00}
\alias{expcub00}
\title{Expectation of CUB distributions}
\usage{
expcub00(m,pai,csi)
}
\arguments{
\item{m}{Number of ordinal categories}

\item{pai}{Uncertainty parameter}

\item{csi}{Feeling parameter}
}
\description{
Compute the expectation of a CUB model without covariates.
}
\examples{
m<-10
pai<-0.3
csi<-0.7
meancub<-expcub00(m,pai,csi)
}
\references{
Piccolo D. (2003). On the moments of a mixture of uniform and shifted binomial random variables.
\emph{Quaderni di Statistica}, \bold{5}, 85--104
}
\seealso{
\code{\link{varcub00}}, \code{\link{expcube}}, \code{\link{varcube}}
}
\keyword{distribution}
