\name{precalc_genlme}
\alias{precalc_genlme}
\title{precalc_genlme}
\usage{
  precalc_genlme(dat, levelnames)
}
\arguments{
  \item{dat}{Dataset with responses, categorize (group),
  covariates and random effects variables. List elements: Y
  is Response variable, F is Categorize (name of group, X
  is Covariates belonging to fixed effects, Z is Covariates
  belonging to random effects.}

  \item{levelnames}{Name of levels in the one-way
  categories.}
}
\value{
  List of precalculated sufficience matrices needed as
  input for the genlme function.
}
\description{
  precalc_genlme calculates sufficient statistics from data
  which is further used in the genlme-function.
}
\details{
  For large amount of data it will be useful to reduce the
  information into sufficience statistics used for making
  inference with Maximum Likelihood estimator.
}
\author{
  Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
}
\references{
  Master Thesis Oyvind Bleka
}
\keyword{LME,EM}

