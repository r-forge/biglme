#' @title genlme
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description genlme estimates parameters a Linear Mixed Effects model which may assume any general covariance structures for the effects.
#' @export
#' @details The Maximum Likelihood are used to estimate the parameters in model using the EM-algorithm. The model assumes gaussian noise with different variance paramater for each group (heterogenous variance structure) and gaussian random effect with any covariance strucuture.
#'
#' The user may specify the restriction for stopping the EM-algorithm using 'eps', minIT and maxIT.
#'
#' The confidence intervals are based on normal-approximated large sample intervals. Details of the EM-algorithm and confidence intervals are found in Masterthesis of Oyvind Bleka.
#'
#' Model: One individual 'j' belongs to a category-level 'i'. Let 'y_ij' be the response, 'X_ij' the covariate-vector for fixed effects, 'Z_ij' the covariate-vector for random effects (equal 1 for random intercept etc.), 'epsilon_ij' is random noise for individual, 'delta_i' is random level effect for category-level 'i'. Then the model is given as 'y_ij=X_ij*beta+Z_ij*delta_i + epsilon_ij'. Here, Cov(epsilon_ij,epsilon_ik)=\{gamma_i for j=k\}\{0 for j!=k\} and Cov(delta_i,delta_l)=[G(phi)^(-1)]_il
#'
#' @param precalc A object returned by function precalc_genlme
#' @param model A list with function-elements' G' and' ldetG' taking a' phi'-parameter where G is the inverse covariance structure of random effects.
#' @param phi0 Startvalues of the' phi'-parameters in the EM-algorithm. 'phi' parameters belonging to the  covariance structure of the random effect levels.
#' @param eps Criterion of euclidean distance for stopping EM-algorithm.
#' @param minIT Minimum number of iterations in the EM-algorithm.
#' @param maxIT Maximum number of iterations in the EM-algorithm.
#' @param hold A list of assumed knowned parameters, given by index numbers. F.ex: list(beta=cbind(1,2),phi=cbind(2,0.1)) means that first beta is fixed=2 and second phi-param is fixed=0.1. 
#' @param theta0 theta0: A list with vector-elements' beta','phi' and' gamma'. Must be same correct size.
#' @return Fitted LME model object
#' \item{est}{Maximum Likelihood estimatation of LME model: 'beta'-Estimated covariate parameters. 'phi'-Estimated parameters to covariance structure of random effects. 'gamma'-Estimated variance parameter to the random noise (one for each categorized level)}
#' \item{OLSest}{Ordinary least squares estimatation of LME model}
#' \item{pred}{E_delta_Y: Mean of random effects. Var_delta_Y: Variance of random effects.}
#' \item{levelnames}{Names of categorized levels.}
#' \item{n_i}{Datasize within each categorized levels.}
#' \item{loglik}{Maximized likelihood value.}
#' \item{iter}{Number of iterations in the EM-algorithm.}
#' \item{timeusage}{Time spent for running function.}
#' \item{modelfit}{AIC and BIC of fitted model.}
#' \item{model}{Same as input}
#' \item{precalc}{Same as input}
#' \item{phi0}{Same as input}
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
#' true=list(beta = rnorm(p,3,1), phi = c(3,1), gamma = rnorm(I,40,1)) #true parameters
#' X = cbind(1,Xsim(p-1,n_i)) #simulate covariate data
#'
#' #Simulate spatial coordinates:
#' XYCRD = cbind(sample(1:I,I),sample(1:I,I)) #specify coordinates for factors
#' XYCRD = XYCRD + matrix(runif(2*I,-0.1,0.1),ncol=2) #add noise
#' rownames(XYCRD) = levelnames
#'
#' #Covariance Prior to level effects:
#' #Specify inverse Covariance matrix as a function of phi:
#' invGam = function(phi) { 
#'  varcov.spatial(XYCRD,cov.model="exponential",func.inv="eigen",cov.pars=c(phi[1],phi[2]),inv=TRUE)$inverse
#' } #invGam(true$phi)
#' #Specify logarithm of determinant of inverse Covariance matrix as a function of phi:
#' logdetG = function(phi) {
#'  -2*varcov.spatial(XYCRD,cov.model="exponential",cov.pars=c(phi[1],phi[2]),det=T)$log.det
#' }  #logdetG(true$phi)
#' modelM2=list(G=invGam,ldetG=logdetG)
#'
#' #Specify random effects:
#' Z = matrix(1,ncol=1,nrow=n) #one intercept effect for each level
#' designM = list(X=X,Z=Z)
#' dat <- genlmesim(model=modelM2,true,designM,n_i,levelnames)
#'
#' #Prefit using simple LME model:
#' lmefitM1 = simplelme(dat,levelnames,phi0=8,eps=10^-5)
#' lmeM1muhat = lmefitM1$pred$E_mu #predicted random effects
#' geoS = as.geodata(cbind(XYCRD,lmeM1muhat)); 
#' mlfit = likfit(geoS,ini.cov.pars=c(6,0.1),fix.nugget=TRUE,cov.model="exponential",messages=FALSE)
#' precalcK1 = precalc_genlme(dat,levelnames) #precalculate for genlme
#' lmefitM2 = genlme(precalcK1,modelM2,phi0=mlfit$cov.pars,eps=10^-5) #fit model
#' lmefitM2CI = genlmeCI(lmefitM2,alpha=0.05,10^-3) #fit (1-alpha)-CI for model
#' }

genlme <-
function(precalc,model,phi0,eps=10^-6,minIT=10^2,maxIT=10^3,hold=NULL,theta0=NULL) {
 #Input: -precalc is a output from function precalc_genlme
 #	 -verbose may be 0: No print, 1:-2loglik, 2:phi in addition, 3:beta in addition
 #       -hold: list of assumed knowned parameters, given by index numbers. 
 #F.ex: list(beta=cbind(1,2),phi=cbind(2,0.1)) means that first beta is fixed=2 and second phi-param is fixed=0.1
 #     -theta0: A list with vector-elements' beta','phi' and' gamma'. Must be same correct size

 #function that gets indices from list2 with values given in list1
 clock = function(time) {
   hours = floor(time[3]/60^2)
   mins = floor( (time[3]-hours*(60^2))/60 )
   secs = round(time[3]%%60)
   return(paste(hours,":",mins,":",secs,"(h:m:s)",sep=""))
 }

 levelnames = precalc$levelnames
 SYTY_i=precalc$SYTY_i
 SXTX_i=precalc$SXTX_i
 SZTZ_i=precalc$SZTZ_i
 SYTX_i=precalc$SYTX_i
 SXTZ_i=precalc$SXTZ_i
 SYTZ_i=precalc$SYTZ_i
 n_i=precalc$n_i
 n = sum(n_i)
 p = dim(SXTX_i[[1]])[2]
 K = dim(SZTZ_i[[1]])[2]
 r = length(phi0)
 I = length(levelnames)
 #Make initialvalues:
 XTX = matrix(0,p,p) #used for LSE
 XTY = matrix(0,nrow=p,ncol=1) #used for LSE
 SSE = rep(0,I) #used for LSE
 for(l in 1:I) {
  XTX = XTX + SXTX_i[[l]]
  XTY = XTY + SYTX_i[l,]
 }
 beta0 = solve(XTX)%*%XTY #use lse as startvalue for beta
 for(l in 1:I) {
  SSE[l] = SYTY_i[l] - 2*SYTX_i[l,]%*%beta0 + t(beta0)%*%SXTX_i[[l]]%*%beta0
 }
 gamma0 = SSE/n_i #using maximum-likelihood-estimatesog sigmasq for each level
 #initialvalues done
 OLSest = c(beta0,rep(NA,r),gamma0)
 if(!is.null(theta0)) {
  beta0 = theta0$beta
  phi0 = theta0$phi
  gamma0 = theta0$gamma
 } 
 beta_k = beta0
 gamma_k = gamma0
 phi_k = phi0
 theta_k = c(beta_k,phi_k,gamma_k) #totally p+r+s parameters
 done = FALSE
 it = 0
 time = system.time({
 while(!done | it<minIT) {
  Gtilde = model$G(phi_k)
  #delta|Y,phi:
  btilde <- rep(NA,(I*K))  
  for(i in 1:I) {
   ind <- c(((1:K)-1)*I+i) #indices of all effect at level i
   Gtilde[ind,ind] <-  Gtilde[ind,ind]  + SZTZ_i[[i]]/gamma_k[i]
   btilde[ind] <- (SYTZ_i[i,]- t(beta_k)%*%SXTZ_i[[i]])/gamma_k[i]
  }
  invGtilde <- solve(Gtilde)
  E_delta_Y = invGtilde%*%btilde
  Var_delta_Y = invGtilde
  EE_delta_Y = Var_delta_Y + E_delta_Y%*%t(E_delta_Y)

  #Beta update:
  WSSX = matrix(0,ncol=p,nrow=p)
  WSXY = matrix(0,ncol=1,nrow=p)
  betaTSXTXTbeta_i = rep(NA,I) #used for calc. likelihood
  SYTXbeta_i = rep(NA,I) #used for calc. likelihood,,,,,,,
  for(i in 1:I) { #go through each level
   betaTSXTXTbeta_i[i] = t(beta_k)%*%SXTX_i[[i]]%*%beta_k
   SYTXbeta_i[i] = SYTX_i[i,]%*%beta_k
   ind <- c(((1:K)-1)*I+i) #indices of all effect at level i
   WSSX = WSSX + SXTX_i[[i]]/gamma_k[i] #.. and sum matrices
   WSXY = WSXY + (SYTX_i[i,]-SXTZ_i[[i]]%*%E_delta_Y[ind])/gamma_k[i]
  }
  beta_kny = as.numeric(solve(WSSX,WSXY)) #new beta given gamma_k
  if(!is.null(hold$beta)) beta_kny[hold$beta[,1]] = hold$beta[,2] #force fixed values  

 #------------------------LIKELIHOOD-----------------------------------#
  loglik = -0.5*( -model$ldetG(phi_k)+n*log(2*pi)+sum(n_i*log(gamma_k))+sum((SYTY_i-2*SYTXbeta_i+betaTSXTXTbeta_i)/gamma_k) + determinant(Gtilde)$mod[1]-t(btilde)%*%invGtilde%*%btilde )
 #---------------------------------------------------------------------#

  #Gamma update: #given updated beta:
  K_i = rep(NA,I)
  for(i in 1:I) {
   ind <- c(((1:K)-1)*I+i) #indices of all effect at level i
   K_i[i] = SYTY_i[i]-2*(SYTX_i[i,]%*%beta_kny+SYTZ_i[i,]%*%E_delta_Y[ind])+ t(beta_kny)%*%SXTX_i[[i]]%*%beta_kny + sum(SZTZ_i[[i]]*EE_delta_Y[ind,ind]) + 2*t(beta_kny)%*%SXTZ_i[[i]]%*%E_delta_Y[ind] 
  }
  gamma_kny =K_i/n_i 
  if(!is.null(hold$gamma)) gamma_kny[hold$gamma[,1]] = hold$gamma[,2] #force fixed values  

  #phi update:
  negQ_phi = function(phi) { 
   phi = exp(phi) #nb: Restriction phi>0 for all elements
   0.5*(-model$ldetG(phi) + sum(model$G(phi)*EE_delta_Y))
  }
  foo = nlm(negQ_phi,log(phi_k))
  phi_kny = exp(foo$est)
  if(!is.null(hold$phi)) phi_kny[hold$phi[,1]] = hold$phi[,2] #force fixed values  

  #verobse:
  d2l = paste(prettyNum(-2*loglik),sep="")
  if(it%%100==0) {
   show(paste("Iteration number ",it,": -2loglik=",d2l,"|phi=",paste(round(phi_k,3),collapse=","),sep=""))
  }
  #if(verbose==1) show(prettyNum(-2*loglik))
  #if(verbose==2) show(c(prettyNum(-2*loglik),prettyNum(phi_kny)))
  #if(verbose==3) show(c(prettyNum(-2*loglik),prettyNum(beta_kny)))
  theta_kny = c(beta_kny,phi_kny,gamma_kny)
  it = it +1
  if((sqrt(sum((theta_kny-theta_k)^2))<eps & it>minIT) | it==maxIT) {
   done=TRUE
  } else {
   beta_k = beta_kny
   phi_k = phi_kny
   gamma_k = gamma_kny
   theta_k = theta_kny
  }  
 } #end iteration k
 })
 show(paste("Timeusage: ",clock(time),". Total iterations=",it, sep=""))
 thetaest = list(beta=beta_kny,phi=phi_kny,gamma=gamma_kny)
 modelfit = list(AIC=-2*loglik+2*length(theta_kny), BIC=-2*loglik+length(theta_kny)*log(n))

 fit = list(est=thetaest,OLSest=OLSest,pred=list(E_delta_Y=E_delta_Y,Var_delta_Y=Var_delta_Y),levelnames=levelnames,n_i=n_i,logLik=loglik, model=model,precalc=precalc,iter=it,phi0=phi0,timeusage=time,modelfit=modelfit)
 return(fit)
}
