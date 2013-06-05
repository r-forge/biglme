\name{genlmesim}
\alias{genlmesim}
\title{genlmesim}
\usage{
  genlmesim(model = NULL, true = NULL, designM = NULL,
    n_i = NULL, levelname = NULL)
}
\arguments{
  \item{model}{A list with function-elements' G' and'
  ldetG' taking a' phi'-parameter.}

  \item{true}{A list with vector-elements' beta','phi' and'
  gamma' used as parameter in the LME model.}

  \item{designM}{A list with covariation matrix' X'
  belonging to fixed effects and random effect matrix' Z'
  belonging to random effects.}

  \item{n_i}{Number of data for each categorizations.}

  \item{levelname}{Name of levels in the one-way
  categories.}
}
\value{
  Dataset with responses, categorize (group), covariates
  and random effects variables \item{Y}{Response variable.}
  \item{F}{Categorize (name of group).} \item{X}{Covariates
  belonging to fixed effects.} \item{Z}{Covariates
  belonging to random effects.} \item{delta}{Sampled random
  effects.}
}
\description{
  genlmesim simulates dataset from a Linear Mixed Effects
  model which may assume any general covariance structures
  for the effects.
}
\details{
  The model assumes gaussian noise with different variance
  paramater for each group (heterogenous variance
  structure) and gaussian random effect with any covariance
  strucuture.
}
\author{
  Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
}
\keyword{LME,EM}

