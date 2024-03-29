% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betar.R
\name{betar}
\alias{betar}
\title{Beta-Binomial distribution}
\usage{
betar(m,csi,phi)
}
\arguments{
\item{m}{Number of ordinal categories}

\item{csi}{Feeling parameter of the Beta-Binomial distribution}

\item{phi}{Overdispersion parameter of the Beta-Binomial distribution}
}
\value{
The vector of length \eqn{m} of the  Beta-Binomial distribution.
}
\description{
Return the Beta-Binomial distribution with parameters \eqn{m}, \eqn{csi} and \eqn{phi}.
}
\examples{
m<-9
csi<-0.8
phi<-0.2
pr<-betar(m,csi,phi)
plot(1:m,pr,type="h", main="Beta-Binomial distribution",xlab="Ordinal categories")
points(1:m,pr,pch=19)
}
\references{
Iannario, M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data,
\emph{Communications in Statistics - Theory and Methods}, \bold{43}, 771--786
}
\seealso{
\code{\link{betabinomial}}
}
\keyword{distribution}
