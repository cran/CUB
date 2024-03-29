% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betabinomialcsi.R
\name{betabinomialcsi}
\alias{betabinomialcsi}
\title{Beta-Binomial probabilities of ordinal responses, given feeling parameter for each observation}
\usage{
betabinomialcsi(m,ordinal,csivett,phi)
}
\arguments{
\item{m}{Number of ordinal categories}

\item{ordinal}{Vector of ordinal responses. Missing values are not allowed: they should be preliminarily deleted
or imputed}

\item{csivett}{Vector of feeling parameters of the Beta-Binomial distribution for given ordinal
responses}

\item{phi}{Overdispersion parameter of the Beta-Binomial distribution}
}
\value{
A vector of the same length as ordinal: each entry is the Beta-Binomial probability for the given observation
for the corresponding feeling and overdispersion parameters.
}
\description{
Compute the Beta-Binomial probabilities of given ordinal responses, with feeling
parameter specified for each observation,
and with the same overdispersion parameter for all the responses.
}
\examples{
data(relgoods)
m<-10
ordinal<-relgoods$Tv
age<-2014-relgoods$BirthYear
no_na<-na.omit(cbind(ordinal,age))
ordinal<-no_na[,1]; age<-no_na[,2]
lage<-log(age)-mean(log(age))
gama<-c(-0.61,-0.31)
phi<-0.16 
csivett<-logis(lage,gama)
pr<-betabinomialcsi(m,ordinal,csivett,phi)
plot(density(pr))
}
\references{
Iannario, M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data,
\emph{Communications in Statistics - Theory and Methods}, \bold{43}, 771--786 \cr
Piccolo D. (2015). Inferential issues for CUBE models with covariates.
\emph{Communications in Statistics - Theory and Methods}, \bold{44}(23), 771--786.
}
\seealso{
\code{\link{betar}}, \code{\link{betabinomial}}
}
\keyword{distribution}
