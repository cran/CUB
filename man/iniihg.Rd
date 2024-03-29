% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iniihg.R
\name{iniihg}
\alias{iniihg}
\title{Moment estimate for the preference parameter of the IHG distribution}
\usage{
iniihg(m,freq)
}
\arguments{
\item{m}{Number of ordinal categories}

\item{freq}{Vector of the absolute frequency distribution of the categories}
}
\value{
Moment estimator of the preference parameter \eqn{\theta}.
}
\description{
Compute the moment estimate of the preference parameter of the IHG distribution.
This preliminary estimate is set as initial value within the optimization procedure for an IHG model
fitting the observed frequencies.
}
\examples{
m<-9
freq<-c(70,51,48,38,29,23,12,10,5)
initheta<-iniihg(m,freq)
}
\references{
D'Elia A. (2003). Modelling ranks using the inverse hypergeometric distribution,
\emph{Statistical Modelling: an International Journal}, \bold{3}, 65--78.
}
\seealso{
\code{\link{inibest}}, \code{\link{inibestcube}}
}
\keyword{htest}
\keyword{utilities}
