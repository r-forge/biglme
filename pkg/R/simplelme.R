#' @title simplelme
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description simplelme makes inference to a one-way Linear Mixed Effects model (assumes indepedent gaussian effect).
#' @export
#' @details The Maximum Likelihood are used to estimate the model using the EM-algorithm. The model assumes gaussian noise with different variance paramater for each group (heterogenous variance structure) and a gaussian independent random intercept effect with constant variance parameter.
#' 
#' The confidence intervals are based on normal-approximated large sample intervals. Details of the EM-algorithm and confidence intervals are found in Masterthesis of Oyvind Bleka.
#'
#' Model: One individual 'j' belongs to a category-level 'i'. Let 'y_ij' be the response, 'X_ij' the covariate-vector for fixed effects, 'mu_i' is random level effect for category-level 'i'. Then the model is given as 'y_ij=X_ij*beta+mu_i + epsilon_ij'. Here, Cov(epsilon_ij,epsilon_ik)={gamma_i for j=k}{0 for j!=k} and Cov(mu_i,mu_l)=\{phi for i=l\}\{0 for i!=l\}.
#'
#' Note that simplelme handels only a special case of LME models which may be fitted using genlme.
#'
#' @param dat Dataset with responses, categorize (group), covariates and random effects variables. List elements: Y is Response variable, F is Categorize (name of group, X is Covariates belonging to fixed effects, Z is Covariates belonging to random effects.
#' @param levelnames  Name of categorized levels.
#' @param phi0 Startvalues of the 'phi'-parameter in the EM-algorithm. Here, 'phi' is the variance parameter of the random effect levels.
#' @param eps Criterion of euclidean distance for stopping EM-algorithm.
#' @param minIT Minimum number of iterations in the EM-algorithm.
#' @param alpha Specifies (1-alpha/2) confidence interval of parameters.
#' @return Fitted simple one-way LME model object
#' \item{est}{Maximum Likelihood estimatation of LME model}
#' \item{OLSest}{Ordinary least squares estimatation of LME model}
#' \item{pred}{E_mu_Y: Mean of random effects. Var_mu_Y: Variance of random effects.}
#' \item{levelnames}{Categorized levelnames.}
#' \item{n_i}{Datasize within each categorized levels.}
#' \item{loglik}{Maximized likelihood value.}
#' \item{iter}{Number of iterations in the EM-algorithm.}
#' \item{timeusage}{Time spent for running function.}
#' \item{modelfit}{AIC and BIC of fitted model.}
#' \item{CI}{Confidence interval of parameters.}
#' @references Master Thesis Oyvind Bleka
#' @keywords LME,EM
#' @examples
#' \dontrun{
#' set.seed(1)
#' require(biglme)
#' require(geoR)
#' Xsim <- function(p,n_i) {
#'  Xtype = sample(1:3,p,replace=TRUE)
#'  Xant = rep(0,3)
#'  for(i in 1:3) Xant[i] = sum(Xtype==i)
#'  X=NULL
#'  I=length(n_i)
#'  n=sum(n_i)
#'  cn_i = c(0,cumsum(n_i))+1 #startindex for each levels
#'  if(p) { #if having any covariables  
#'   for(i in 1:I) { #for each levels
#'    Xi = matrix(NA,nrow=n_i[i],ncol=p) 
#'    Xi[,which(Xtype==1)] = matrix(rnorm(n_i[i]*Xant[1],5,3),nrow=n_i[i],ncol=Xant[1])
#'   Xi[,which(Xtype==2)] = matrix(rpois(n_i[i]*Xant[2],3),nrow=n_i[i],ncol=Xant[2])
#'    Xi[,which(Xtype==3)] = matrix(rbinom(n_i[i]*Xant[3],1,0.3),nrow=n_i[i],ncol=Xant[3])
#'    if(i==1) { X = Xi 
#'    } else { X = rbind(X,Xi) } #just add up the matrix
#'   }
#'  }
#'  return(X)
#' }
#' I = 30 #number of effectlevels: 
#' levelnames = paste("place",1:I,sep="") #name of levels
#' nlvl = 1000 #expected number of data per level
#' n=I*nlvl #total number of data
#' n_i = c(rmultinom(1,n,runif(I,0.3,0.7))) #gen number of observations at each level
#' p = 4 #number of covariates
#' true=list(beta = rnorm(p,3,1), phi = c(3), gamma = rnorm(I,40,1)) #true parameters
#' X = cbind(1,Xsim(p-1,n_i)) #simulate covariate data
#'
#' #Covariance Prior to level effects:
#' invGam = function(phi) { 
#'  diag(rep(1,I))/phi
#' } #invGam(true$phi)
#' #Specify logarithm of determinant of inverse Covariance matrix as a function of phi:
#' logdetG = function(phi) {
#'  -I*phi
#' } 
#' modelM1=list(G=invGam,ldetG=logdetG)
#' Z = matrix(1,ncol=1,nrow=n) #one intercept effect for each level
#' designM = list(X=X,Z=Z)
#'
#' dat <- genlmesim(model=modelM1,true,designM,n_i,levelnames)
#' lmefitM1 = simplelme(dat,levelnames,phi0=8,eps=10^-5)
#' }

simplelme <-
function(dat,levelnames,phi0=1,eps=10^-5,minIT=10^2,alpha=0.05) {

 #function that gets indices from list2 with values given in list1
 compareIndices = function(list1,list2) { 
  index = rep(NA,length(list1))
  for(i in 1:length(list1)) {
   index[i] = which(list1[i]==list2)
  }
  return(index)
 }
 clock = function(time) {
   hours = floor(time[3]/60^2)
   mins = floor( (time[3]-hours*(60^2))/60 )
   secs = round(time[3]%%60)
   return(paste(hours,":",mins,":",secs,"(h:m:s)",sep=""))
 }

 time = system.time( {
 p = dim(dat$X)[2]
 I = length(levelnames)
 #precalculations:
 sumsq = function(x){ sum(x^2)} #returns sumsquares
 agg = aggregate(dat$Y,by=list(dat$F),length)
 agglvls = agg[,1]
 gind = compareIndices(levelnames,agglvls) #make agg.ed data in right order.
 #agglvls[gind]==levelnames #OK. Now use index on all agg.ed-data
 n_i = agg[gind,2]
 S_i = aggregate(dat$Y,by=list(dat$F),sum)[gind,2]
 SS_i = aggregate(dat$Y,by=list(dat$F),sumsq)[gind,2]
 SYX_i = as.matrix(aggregate(dat$Y*dat$X,by=list(dat$F),sum)[gind,2:(1+p)]) #(Ixp) matrix
 SX_i = as.matrix(aggregate(dat$X,by=list(dat$F),sum)[gind,2:(1+p)]) #(Ixp) matrix
 SXTX_i = matrix(list(),ncol=1,nrow=I) #list with (pxp) matrices 
 XTX = matrix(0,p,p) #used for LSE
 XTY = matrix(0,nrow=p,ncol=1) #used for LSE
 SSE = rep(0,I) #used for LSE
 for(l in 1:I) {
  S = matrix(0,ncol=p,nrow=p)
  subX = matrix(dat$X[dat$F==levelnames[l],],ncol=p) #take out subset of current level
  SXTX_i[[l]] = t(subX)%*%subX
  XTX = XTX + SXTX_i[[l]]
  XTY = XTY + SYX_i[l,]
 }
 beta0 = solve(XTX)%*%XTY #use lse as startvalue for beta
 for(l in 1:I) {
  SSE[l] = SS_i[l] - 2*SYX_i[l,]%*%beta0 + t(beta0)%*%SXTX_i[[l]]%*%beta0
 }
 gamma0 = SSE/n_i #using maximum-likelihood-estimates 
 #END precalculations#
 beta_k = beta0
 gamma_k = gamma0
 phi_k = phi0
 theta_k = c(beta_k,phi_k,gamma_k) #totally p+r+s parameters
 
 done = FALSE
 it = 0
 while(!done | it<minIT) {
  E_mu_Y = as.numeric(phi_k*(S_i-SX_i%*%beta_k)/(gamma_k+n_i*phi_k)) #E[mu_i|Y,th^k]
  Var_mu_Y = phi_k*gamma_k/(gamma_k+n_i*phi_k) #Var[mu_i|Y,th^k]
  E2_mu_Y = Var_mu_Y + E_mu_Y^2 #E[mu_i^2|Y,th^k]

  #update parameters: 
  WSSX = matrix(0,ncol=p,nrow=p)
  for(i in 1:I) { #go through each level
    WSSX = WSSX + SXTX_i[[i]]/gamma_k[i] #.. and sum matrices
  }
  beta_kny = as.numeric(solve(WSSX)%*%colSums(diag(c(1/gamma_k))%*%(SYX_i-SX_i*E_mu_Y))) 

  #precalcs. for gamma:
  sumsqXB = rep(NA,I)
  for(i in 1:I) {
   sumsqXB[i] = t(beta_kny)%*%SXTX_i[[i]]%*%beta_kny
  }
  gamma_kny = (SS_i-2*c(SYX_i%*%beta_kny)-2*S_i*E_mu_Y+sumsqXB+2*E_mu_Y*c(SX_i%*%beta_kny) +n_i*E2_mu_Y)/n_i  

  phi_kny = sum(E2_mu_Y)/I

  it = it +1
  theta_kny = c(beta_kny,phi_kny,gamma_kny)
  err = sqrt(sum((theta_k-theta_kny)^2))
  if(err<eps[1]) {
   done=TRUE
  } else {
   beta_k = beta_kny
   phi_k = phi_kny
   gamma_k = gamma_kny
   theta_k = theta_kny
  }
  if(it%%10^3==0) show(paste("Iterations: ",it,sep=""))
 }

 loglik = function(theta) {
  beta=theta[1:p]
  phi=theta[p+1]
  gamma=theta[(p+2):(p+1+I)]
  lambda = gamma+n_i*phi  
  BXTXB = rep(NA,I)
  for(i in 1:I) {
   BXTXB[i] = t(beta)%*%SXTX_i[[i]]%*%beta
  }
  genres = SS_i-2*c(SYX_i%*%beta)+BXTXB-phi*lambda^(-1)*(S_i-c(SX_i%*%beta))^2
  lvlval = n_i*log(2*pi) + (n_i-1)*log(gamma) + log(lambda) + gamma^(-1)*genres
  loglik = -0.5*sum(lvlval)
  return(loglik)
 }

 #Calculate confidence interval:
 #Calc. observed information matrix: 
 hessian = function(theta) {
  beta=theta[1:p]
  phi=theta[p+1]
  gamma=theta[(p+2):(p+1+I)]
  lambda = (gamma+n_i*phi) #substitute
  l_beta2 = matrix(0,p,p)
  BXTXB = rep(NA,I)
  for(i in 1:I) {
   l_beta2=l_beta2 + gamma[i]^(-1)*(-SXTX_i[[i]]+phi*lambda[i]^(-1)*SX_i[i,]%*%t(SX_i[i,]))
   BXTXB[i] = t(beta)%*%SXTX_i[[i]]%*%beta
  }
  g_i = SS_i - 2*SYX_i%*%beta + BXTXB - phi*lambda^(-1)*(S_i-SX_i%*%beta)^2
  l_phi2 = 0.5*sum( lambda^(-2)*n_i^2 - 2*gamma^(-1)*(S_i-SX_i%*%beta)^2*( lambda^(-2)*n_i -phi*lambda^(-3)*n_i^2 ) )
  l_gamma2 = -0.5*(-(n_i-1)*gamma^(-2) - lambda^(-2) + 2*gamma^(-3)*g_i - 2*phi*(S_i-SX_i%*%beta)^2*(gamma^(-2)*lambda^(-2) + gamma^(-1)*lambda^(-3) ) )
  l_betagamma = matrix(0,nrow=I,ncol=p)
  l_betaphi = matrix(0,nrow=I,ncol=p)
  l_phigamma = rep(NA,I)
  for(i in 1:I) {
   f_i = S_i[i]*SX_i[i,] - t(SX_i[i,])*c(SX_i[i,]%*%beta)
   l_betagamma[i,] = -gamma[i]^(-2)*( SYX_i[i,]-beta%*%SXTX_i[[i]] - phi*lambda[i]^(-1)*f_i) + gamma[i]^(-1)*phi*lambda[i]^(-2)*f_i 
   l_betaphi[i,] = gamma[i]^(-1)*f_i*( lambda[i]^(-1)-phi*lambda[i]^(-2)*n_i[i] )
  }
  l_betaphi = -colSums(l_betaphi)
  l_phigamma = 0.5*(lambda^(-2)*n_i - (S_i-SX_i%*%beta)^2*( gamma^(-2)*(lambda^(-1)-phi*lambda^(-2)*n_i)+gamma^(-1)*(lambda^(-2)-2*phi*lambda^(-3)*n_i) ) )
 
  H = matrix(0,(p+1+I),(p+1+I)) #must be symmetrical!
  H[1:p,1:p] = l_beta2
  H[p+1,p+1] = l_phi2
  diag(H[(p+2):(p+1+I),(p+2):(p+1+I)]) = l_gamma2
  H[p+1,1:p] <- H[1:p,p+1] <- l_betaphi
  H[(p+2):(p+1+I),1:p] <- l_betagamma
  H[1:p,(p+2):(p+1+I)] <- t(l_betagamma)
  H[p+1,(p+2):(p+1+I)] <- l_phigamma
  H[(p+2):(p+1+I),p+1] <- t(l_phigamma)
  return(H)
 }
 loglikval = loglik(theta_kny) 
 modelfit = list(AIC=-2*loglikval+2*length(theta_kny), BIC=-2*loglikval+length(theta_kny)*log(sum(n_i)))

 J = hessian(theta_kny) #observed information matrix
 stderr = sqrt(diag(solve(-J)))
 CI_N = cbind(theta_kny + qnorm(alpha/2)*stderr,theta_kny,theta_kny + qnorm(1-alpha/2)*stderr)
 colnames(CI_N) = c(paste("(",alpha/2,"%",sep=""),"est",paste(1-alpha/2,"%)",sep=""))
 rownames(CI_N)=c(paste("beta",1:p,sep=""),paste("phi",1,sep=""),paste("gamma",1:length(levelnames),sep=""))
 })
 show(paste("Timeusage: ",clock(time)," with number of iterations=",it, sep=""))

 MLest = list()
 MLest$beta = theta_kny[1:length(beta0)]
 MLest$phi = theta_kny[length(beta0)+1]
 MLest$gamma = theta_kny[(length(beta0)+2):length(theta_kny)]
 fit = list(est=MLest,OLSest=c(beta0,NA,gamma0),pred=list(E_mu_Y=E_mu_Y,Var_mu_Y=Var_mu_Y),logLik=loglikval,hessian=J,CI=CI_N,iter=it,levelnames=levelnames, n_i=n_i,timeusage=time,modelfit=modelfit)
 return(fit)
}
