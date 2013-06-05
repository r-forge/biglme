\name{genlme}
\alias{genlme}
\title{genlme}
\usage{
  genlme(precalc, model, phi0, eps = 10^-6, minIT = 10^2,
    maxIT = 10^3, hold = NULL, theta0 = NULL)
}
\arguments{
  \item{precalc}{A object returned by function
  precalc_genlme}

  \item{model}{A list with function-elements' G' and'
  ldetG' taking a' phi'-parameter where G is the inverse
  covariance structure of random effects.}

  \item{phi0}{Startvalues of the' phi'-parameters in the
  EM-algorithm. 'phi' parameters belonging to the
  covariance structure of the random effect levels.}

  \item{eps}{Criterion of euclidean distance for stopping
  EM-algorithm.}

  \item{minIT}{Minimum number of iterations in the
  EM-algorithm.}

  \item{maxIT}{Maximum number of iterations in the
  EM-algorithm.}

  \item{hold}{A list of assumed knowned parameters, given
  by index numbers. F.ex:
  list(beta=cbind(1,2),phi=cbind(2,0.1)) means that first
  beta is fixed=2 and second phi-param is fixed=0.1.}

  \item{theta0}{theta0: A list with vector-elements'
  beta','phi' and' gamma'. Must be same correct size.}
}
\value{
  Fitted LME model object \item{est}{Maximum Likelihood
  estimatation of LME model: 'beta'-Estimated covariate
  parameters. 'phi'-Estimated parameters to covariance
  structure of random effects. 'gamma'-Estimated variance
  parameter to the random noise (one for each categorized
  level)} \item{OLSest}{Ordinary least squares estimatation
  of LME model} \item{pred}{E_delta_Y: Mean of random
  effects. Var_delta_Y: Variance of random effects.}
  \item{levelnames}{Names of categorized levels.}
  \item{n_i}{Datasize within each categorized levels.}
  \item{loglik}{Maximized likelihood value.}
  \item{iter}{Number of iterations in the EM-algorithm.}
  \item{timeusage}{Time spent for running function.}
  \item{modelfit}{AIC and BIC of fitted model.}
  \item{model}{Same as input} \item{precalc}{Same as input}
  \item{phi0}{Same as input}
}
\description{
  genlme estimates parameters a Linear Mixed Effects model
  which may assume any general covariance structures for
  the effects.
}
\details{
  The Maximum Likelihood are used to estimate the
  parameters in model using the EM-algorithm. The model
  assumes gaussian noise with different variance paramater
  for each group (heterogenous variance structure) and
  gaussian random effect with any covariance strucuture.

  The user may specify the restriction for stopping the
  EM-algorithm using 'eps', minIT and maxIT.

  The confidence intervals are based on normal-approximated
  large sample intervals. Details of the EM-algorithm and
  confidence intervals are found in Masterthesis of Oyvind
  Bleka.

  Model: One individual 'j' belongs to a category-level
  'i'. Let 'y_ij' be the response, 'X_ij' the
  covariate-vector for fixed effects, 'Z_ij' the
  covariate-vector for random effects (equal 1 for random
  intercept etc.), 'epsilon_ij' is random noise for
  individual, 'delta_i' is random level effect for
  category-level 'i'. Then the model is given as
  'y_ij=X_ij*beta+Z_ij*delta_i + epsilon_ij'. Here,
  Cov(epsilon_ij,epsilon_ik)=\{gamma_i for j=k\}\{0 for
  j!=k\} and Cov(delta_i,delta_l)=[G(phi)^(-1)]_il
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
true=list(beta = rnorm(p,3,1), phi = c(3,1), gamma = rnorm(I,40,1)) #true parameters
X = cbind(1,Xsim(p-1,n_i)) #simulate covariate data

#Simulate spatial coordinates:
XYCRD = cbind(sample(1:I,I),sample(1:I,I)) #specify coordinates for factors
XYCRD = XYCRD + matrix(runif(2*I,-0.1,0.1),ncol=2) #add noise
rownames(XYCRD) = levelnames

#Covariance Prior to level effects:
#Specify inverse Covariance matrix as a function of phi:
invGam = function(phi) {
 varcov.spatial(XYCRD,cov.model="exponential",func.inv="eigen",cov.pars=c(phi[1],phi[2]),inv=TRUE)$inverse
} #invGam(true$phi)
#Specify logarithm of determinant of inverse Covariance matrix as a function of phi:
logdetG = function(phi) {
 -2*varcov.spatial(XYCRD,cov.model="exponential",cov.pars=c(phi[1],phi[2]),det=T)$log.det
}  #logdetG(true$phi)
modelM2=list(G=invGam,ldetG=logdetG)

#Specify random effects:
Z = matrix(1,ncol=1,nrow=n) #one intercept effect for each level
designM = list(X=X,Z=Z)
dat <- genlmesim(model=modelM2,true,designM,n_i,levelnames)

#Prefit using simple LME model:
lmefitM1 = simplelme(dat,levelnames,phi0=8,eps=10^-5)
lmeM1muhat = lmefitM1$pred$E_mu #predicted random effects
geoS = as.geodata(cbind(XYCRD,lmeM1muhat));
mlfit = likfit(geoS,ini.cov.pars=c(6,0.1),fix.nugget=TRUE,cov.model="exponential",messages=FALSE)
precalcK1 = precalc_genlme(dat,levelnames) #precalculate for genlme
lmefitM2 = genlme(precalcK1,modelM2,phi0=mlfit$cov.pars,eps=10^-5) #fit model
lmefitM2CI = genlmeCI(lmefitM2,alpha=0.05,10^-3) #fit (1-alpha)-CI for model
}
}
\author{
  Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
}
\references{
  Master Thesis Oyvind Bleka
}
\keyword{LME,EM}

