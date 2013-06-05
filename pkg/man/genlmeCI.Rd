\name{genlmeCI}
\alias{genlmeCI}
\title{genlmeCI}
\usage{
  genlmeCI(lmefit, alpha = 0.05, eps = 10^-5)
}
\arguments{
  \item{lmefit}{A fitted object returned from function
  genlme.}

  \item{alpha}{Specifies (1-alpha/2) confidence interval of
  parameters.}

  \item{eps}{Used for numerical first and second partial
  derivation of change of the inverse covariance function
  on the phi-parameter.}
}
\value{
  A List with the following \item{CI}{Confidence interval }
  \item{hessian}{The estimated Hessian matrix to the
  parameters.} \item{Var_thetahat}{Large-sample variance of
  the parameterestimators.} \item{timeusage}{Time spent for
  calculating confidence interval}
}
\description{
  genlmeCI calculates confidence interval of parameters in
  a Linear Mixed Effects model which may assume any general
  covariance structures for the effects.
}
\details{
  genlmeCI takes a Linear Mixed Effect fitted object
  returned from function genlme.
}
\author{
  Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
}
\references{
  Master Thesis Oyvind Bleka
}
\keyword{Confidence}
\keyword{interval.}
\keyword{LME,EM,Large}
\keyword{Sample}

