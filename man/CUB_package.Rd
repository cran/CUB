% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CUB_package.R
\name{CUB_package}
\alias{CUB_package}
\title{CUB package}
\description{
The analysis of human perceptions is often carried out by resorting to questionnaires,
where respondents are asked to express ratings about the items being evaluated. The standard goal of the
statistical framework proposed for this kind of data (e.g. cumulative models) is to explicitly characterize
the respondents' perceptions about a latent trait, by taking into account, at the same time,
the ordinal categorical scale of measurement of the   involved statistical variables.\cr
The new class of models starts from a particular assumption about the unconscious mechanism leading individuals' responses
to choose an ordinal category on a rating scale. The basic idea derives from the awareness that two latent
components move the psychological process of selection among discrete alternatives: attractiveness
towards the item and uncertainty in the response. Both components of models concern the stochastic
mechanism in term of feeling, which is an internal/personal movement of the subject towards the item,
and uncertainty pertaining to the final choice.\cr
Thus, on the basis of experimental data and statistical motivations, the response distribution is modelled
as the convex Combination of a discrete Uniform and a shifted Binomial random variable (denoted as CUB model)
whose parameters may be consistently estimated and validated by maximum likelihood inference.
In addition, subjects' and objects' covariates can be included in the model in order to assess how the
characteristics of the respondents may affect the ordinal score. \cr
CUB models have been firstly introduced by Piccolo (2003) and implemented on real datasets concerning ratings and rankings
by D'Elia and Piccolo (2005). \cr
The CUB package allows the user to estimate and test CUB models and their extensions by using maximum
likelihood methods: see Piccolo and Simone (2019a, 2019b) for an updated overview of methodological developments and applications.
The accompanying vignettes supplies the user with detailed usage instructions and examples. \cr
Acknowledgements: The Authors are grateful to Maria Antonietta Del Ferraro, Francesco Miranda and
Giuseppe Porpora for their preliminary support in the implementation of the first version of the package.
}
\details{
\tabular{ll}{
Package: \tab CUB\cr
Type: \tab Package\cr
Version: \tab 1.1.4\cr
Date: \tab 2017-10-11\cr
License: GPL-2 | GPL-3
}
}
\references{
D'Elia A. (2003). Modelling ranks using the inverse hypergeometric distribution,
\emph{Statistical Modelling: an International Journal}, \bold{3}, 65--78 \cr
Piccolo D. (2003). On the moments of a mixture of uniform and shifted binomial random variables,
\emph{Quaderni di Statistica}, \bold{5}, 85--104 \cr
D'Elia A. and Piccolo D. (2005).  A mixture model for preferences data analysis,
\emph{Computational Statistics & Data Analysis},  \bold{49}, 917--937 \cr
Piccolo D. and Simone R. (2019a).  The class of CUB models: statistical foundations, inferential issues and empirical evidence.
\emph{Statistical Methods and Applications}, \bold{28}(3), 389--435.\cr
Piccolo D. and Simone R. (2019b). Rejoinder to the discussions: The class of CUB models: statistical foundations, inferential issues and empirical evidence.
\emph{Statistical Methods and Applications}, \bold{28}(3), 477-493.\cr
Capecchi S. and Piccolo D. (2017). Dealing with heterogeneity in ordinal responses,
\emph{Quality and Quantity}, \bold{51}(5), 2375--2393 \cr
\emph{Metron}, \bold{74}(2), 233--252.\cr
Iannario M. and Piccolo D. (2016b). A generalized framework for modelling ordinal data.
\emph{Statistical Methods and Applications}, \bold{25}, 163--189.\cr
}
\author{
Maria Iannario, Domenico Piccolo, Rosaria Simone
}
\keyword{package}
