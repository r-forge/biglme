\name{simplelme}
\alias{simplelme}
\title{simplelme}
\usage{
  simplelme(dat, levelnames, phi0 = 1, eps = 10^-5,
    minIT = 10^2, alpha = 0.05)
}
\arguments{
  \item{dat}{Dataset with responses, categorize (group),
  covariates and random effects variables. List elements: Y
  is Response variable, F is Categorize (name of group, X
  is Covariates belonging to fixed effects, Z is Covariates
  belonging to random effects.}

  \item{levelnames}{Name of categorized levels.}

  \item{phi0}{Startvalues of the 'phi'-parameter in the
  EM-algorithm. Here, 'phi' is the variance parameter of
  the random effect levels.}

  \item{eps}{Criterion of euclidean distance for stopping
  EM-algorithm.}

  \item{minIT}{Minimum number of iterations in the
  EM-algorithm.}

  \item{alpha}{Specifies (1-alpha/2) confidence interval of
  parameters.}
}
\value{
  Fitted simple one-way LME model object \item{est}{Maximum
  Likelihood estimatation of LME model}
  \item{OLSest}{Ordinary least squares estimatation of LME
  model} \item{pred}{E_mu_Y: Mean of random effects.
  Var_mu_Y: Variance of random effects.}
  \item{levelnames}{Categorized levelnames.}
  \item{n_i}{Datasize within each categorized levels.}
  \item{loglik}{Maximized likelihood value.}
  \item{iter}{Number of iterations in the EM-algorithm.}
  \item{timeusage}{Time spent for running function.}
  \item{modelfit}{AIC and BIC of fitted model.}
  \item{CI}{Confidence interval of parameters.}
}
\description{
  simplelme makes inference to a one-way Linear Mixed
  Effects model (assumes indepedent gaussian effect).
}
\details{
  The Maximum Likelihood are used to estimate the model
  using the EM-algorithm. The model assumes gaussian noise
  with different variance paramater for each group
  (heterogenous variance structure) and a gaussian
  independent random intercept effect with constant
  variance parameter.

  The confidence intervals are based on normal-approximated
  large sample intervals. Details of the EM-algorithm and
  confidence intervals are found in Masterthesis of Oyvind
  Bleka.

  Model: One individual 'j' belongs to a category-level
  'i'. Let 'y_ij' be the response, 'X_ij' the
  covariate-vector for fixed effects, 'mu_i' is random
  level effect for category-level 'i'. Then the model is
  given as 'y_ij=X_ij*beta+mu_i + epsilon_ij'. Here,
  Cov(epsilon_ij,epsilon_ik)={gamma_i for j=k}{0 for j!=k}
  and Cov(mu_i,mu_l)=\{phi for i=l\}\{0 for i!=l\}.

  Note that simplelme handels only a special case of LME
  models which may be fitted using genlme.
}
\examples{
\dontrun{
set.seed(1)
require(biglme)
require(geoR)
Xsim <- function(p,n_i) {
 Xtype = sample(1:3,p,replace=TRUE)
 Xant = rep(0,3)
 for(i in 1:3) Xant[i] = sum(Xtype==i)
 X=NULL
 I=length(n_i)
 n=sum(n_i)
 cn_i = c(0,cumsum(n_i))+1 #startindex for each levels
 if(p) { #if having any covariables
  for(i in 1:I) { #for each levels
   Xi = matrix(NA,nrow=n_i[i],ncol=p)
   Xi[,which(Xtype==1)] = matrix(rnorm(n_i[i]*Xant[1],5,3),nrow=n_i[i],ncol=Xant[1])
  Xi[,which(Xtype==2)] = matrix(rpois(n_i[i]*Xant[2],3),nrow=n_i[i],ncol=Xant[2])
   Xi[,which(Xtype==3)] = matrix(rbinom(n_i[i]*Xant[3],1,0.3),nrow=n_i[i],ncol=Xant[3])
   if(i==1) { X = Xi
   } else { X = rbind(X,Xi) } #just add up the matrix
  }
 }
 return(X)
}
I = 30 #number of effectlevels:
levelnames = paste("place",1:I,sep="") #name of levels
nlvl = 1000 #expected number of data per level
n=I*nlvl #total number of data
n_i = c(rmultinom(1,n,runif(I,0.3,0.7))) #gen number of observations at each level
p = 4 #number of covariates
true=list(beta = rnorm(p,3,1), phi = c(3), gamma = rnorm(I,40,1)) #true parameters
X = cbind(1,Xsim(p-1,n_i)) #simulate covariate data

#Covariance Prior to level effects:
invGam = function(phi) {
 diag(rep(1,I))/phi
} #invGam(true$phi)
#Specify logarithm of determinant of inverse Covariance matrix as a function of phi:
logdetG = function(phi) {
 -I*phi
}
modelM1=list(G=invGam,ldetG=logdetG)
Z = matrix(1,ncol=1,nrow=n) #one intercept effect for each level
designM = list(X=X,Z=Z)

dat <- genlmesim(model=modelM1,true,designM,n_i,levelnames)
lmefitM1 = simplelme(dat,levelnames,phi0=8,eps=10^-5)
}
}
\author{
  Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
}
\references{
  Master Thesis Oyvind Bleka
}
\keyword{LME,EM}

