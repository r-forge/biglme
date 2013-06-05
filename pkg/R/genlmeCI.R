#' @title genlmeCI
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description genlmeCI calculates confidence interval of parameters in a Linear Mixed Effects model which may assume any general covariance structures for the effects.
#' @export
#' @details genlmeCI takes a Linear Mixed Effect fitted object returned from function genlme.
#' 
#' @param lmefit A fitted object returned from function genlme.
#' @param alpha Specifies (1-alpha/2) confidence interval of parameters.
#' @param eps Used for numerical first and second partial derivation of change of the inverse covariance function on the phi-parameter. 
#' @return A List with the following
#' \item{CI}{Confidence interval }
#' \item{hessian}{The estimated Hessian matrix to the parameters.}
#' \item{Var_thetahat}{Large-sample variance of the parameterestimators.}
#' \item{timeusage}{Time spent for calculating confidence interval}
#' @references Master Thesis Oyvind Bleka
#' @keywords LME,EM,Large Sample Confidence interval.

genlmeCI <-
function(lmefit,alpha=0.05,eps=10^-5) {
 #lmefit is a fitted object from lmeMgen-routine
 #precalc is precalculated values from precalc_lmeMgen
 #calc. CI using likelihood-derivatives of est. theta from lmefit.
 #If thetapar is specified, only the hessian is computed
 #Function that gives numerical derivatives and hessian matrix
 #of a given G-function
 deltaGlist = function(G,phi,eps) {
  r = length(phi)
  f0 = G(phi)
  DG = rep(list(),r)
  DDG = matrix(list(),r,r)
  fval = matrix(list(),nrow=r,ncol=2) #saves temp. values
  ffval = matrix(list(),nrow=r,ncol=r) #saves temp. values
  zero = rep(0,r)
  #first derivatives:
  for(k in 1:r) {
   e = zero
   e[k] = eps
   fval[[k,1]] = G(phi+e)
   fval[[k,2]] = G(phi-e)
   DG[[k]] = (fval[[k,1]]-f0)/eps   #insert derivative
  }
  #second derivatives:
  for(k in 1:r) {
   e = zero
   e[k] = eps
   for(l in k:r) {
    if(k!=l) {
     ee = e #temp
     ee[l] = eps #place in addition
     ffval = G(phi+ee)
     #make symmetri for cross-derivative:
     DDG[[k,l]] <- DDG[[l,k]] <- (ffval-fval[[k,1]]-fval[[l,1]]+f0)/eps^2
    } else {
     DDG[[k,k]] = (fval[[k,1]]-2*f0+fval[[k,2]])/eps^2
    }
   }
  }
  return(list(deltaG=DG,deltasqG=DDG))
 }

 clock = function(time) {
   hours = floor(time[3]/60^2)
   mins = floor( (time[3]-hours*(60^2))/60 )
   secs = round(time[3]%%60)
   return(paste(hours,":",mins,":",secs,"(h:m:s)",sep=""))
 }

 #Function starts
 precalc = lmefit$precalc
 E_delta_Y = lmefit$pred$E_delta_Y
 Var_delta_Y = round(lmefit$pred$Var_delta_Y,15)
 EE_delta_Y = Var_delta_Y + E_delta_Y%*%t(E_delta_Y)
 beta = lmefit$est$beta
 phi = lmefit$est$phi
 gamma = lmefit$est$gamma
 theta=c(beta,phi,gamma)
 p = length(beta)
 r = length(phi)
 n_i=precalc$n_i #number of obs. for each levels
 levelnames = lmefit$levelnames
 I = length(levelnames) 
 SYTY_i=precalc$SYTY_i
 SXTX_i=precalc$SXTX_i
 SZTZ_i=precalc$SZTZ_i
 SYTX_i=precalc$SYTX_i
 SXTZ_i=precalc$SXTZ_i
 SYTZ_i=precalc$SYTZ_i
 K = dim(SZTZ_i[[1]])[2] #number of effects

 #dependency function:
 G = lmefit$model$G(phi)
 Ginv = solve(G)
 dGlist = deltaGlist(lmefit$model$G,phi,eps) #get the derivatives-matrices
 deltaG = dGlist$deltaG #r long list with derivatives
 deltasqG = dGlist$deltasqG #rxr -matrixlist with second derivatives

 time = system.time( {
 #precalc for E_deltasq-matrix:
  WSSX = matrix(0,ncol=p,nrow=p)
  K_i = rep(NA,I)
  for(i in 1:I) { #go through each level
   ind <- c(((1:K)-1)*I+i) #indices of all effect at level i
   WSSX = WSSX + SXTX_i[[i]]/gamma[i] #.. and sum matrices
   K_i[i] = SYTY_i[i]-2*(SYTX_i[i,]%*%beta+SYTZ_i[i,]%*%E_delta_Y[ind])+ t(beta)%*%SXTX_i[[i]]%*%beta + sum(SZTZ_i[[i]]*EE_delta_Y[ind,ind]) + 2*t(beta)%*%SXTZ_i[[i]]%*%E_delta_Y[ind] 
  }
 #Expectation of second derivatives of complete likelihood#
  Edeltasql_beta2 = -WSSX 
  Edeltasql_gamma2 = -0.5*( -n_i*gamma^(-2)+2*gamma^(-3)*K_i )
  Edeltasql_betagamma = matrix(0,ncol=I,nrow=p)
  for(i in 1:I) { 
   ind <- c(((1:K)-1)*I+i) #indices of all effect at level i
   Edeltasql_betagamma[,i] = -gamma[i]^(-2)*(SYTX_i[i,]-c(SXTZ_i[[i]]%*%E_delta_Y[ind])-t(beta)%*%SXTX_i[[i]])
  }
  Edeltasql_phi2 = matrix(0,r,r) #hessian for phi
  for(i in 1:r) 
  for(j in 1:r) {
   deltaGi = deltaG[[i]]
   deltaGj = deltaG[[j]]
   deltasqGij = deltasqG[[i,j]]
   Edeltasql_phi2[i,j] = -0.5*(sum((Ginv%*%deltaGj)*(Ginv%*%deltaGi))-sum(Ginv*deltasqGij)+sum(deltasqGij*EE_delta_Y) )
  }
  Edeltasql = rbind( 
   cbind(Edeltasql_beta2, matrix(0,nrow=p,ncol=r) , Edeltasql_betagamma),
   cbind(matrix(0,nrow=r,ncol=p),Edeltasql_phi2,matrix(0,nrow=r,ncol=I)),
   cbind(t(Edeltasql_betagamma),matrix(0,nrow=I,ncol=r),diag(Edeltasql_gamma2))
  )
 #END: Expectation of second derivatives#
  
 #Analytical Var_score:
  #Precalcs:
  #xi: rearranging delta to be sorted after levels.
  #let xi~N(m,Psi) <=> rearranged delta|Y ~N(E_delta_Y,Var_delta_Y) after levels
  COV_ddTZTZd = matrix(NA,nrow=(I*K),ncol=I) #cubic: COV[delta,deltajTZjTZjdeltaj] qx1 for each level
  COV_ddTdeltaGld = matrix(NA,nrow=(I*K),ncol=r) #cubic: COV[delta,deltaTDGphiDkdelta] qxr
  COV_dTZTZddTdeltaGld = matrix(0,nrow=I,ncol=r) #COV[deltajTZjTZjdeltaj,deltaTDGphiDkdelta] Ixr
  COV_dTZTZddTZTZd = matrix(0,I,I) #COV[deltaiTZiTZideltai,deltajTZjTZjdeltaj] IxI

  #dTZTZddTZTZd
  for(i in 1:I) {
   indi <- c(((1:K)-1)*I+i) #indices of all effect at level i
   for(j in i:I) {
    indj <- c(((1:K)-1)*I+j) #indices of all effect at level j
    Trase = sum(diag(SZTZ_i[[i]]%*%Var_delta_Y[indi,indj]%*%SZTZ_i[[j]]%*%Var_delta_Y[indj,indi]))
    COV_dTZTZddTZTZd[i,j] = 2*Trase + 4*t(E_delta_Y[indi])%*%SZTZ_i[[i]]%*%Var_delta_Y[indi,indj]%*%SZTZ_i[[j]]%*%E_delta_Y[indj]
   }
  } # COV_dTZTZddTZTZd should be symmetric: There are some round off-errors

  #fill in symmetry:
  for(i in 1:(I-1)) {
   for(j in (i+1):I) {
    COV_dTZTZddTZTZd[j,i] = COV_dTZTZddTZTZd[i,j]
   }
  }

  #COV_ddTZTZd:
  for(i in 1:I) {
   indi <- c(((1:K)-1)*I+i) #indices of all effect at level i
   TEMP = SZTZ_i[[i]]%*%E_delta_Y[indi] #Kx1 - matrix
   for(j in 1:I) {
    indj <- c(((1:K)-1)*I+j) #indices of all effect at level j
    COV_ddTZTZd[indj,i] = 2*Var_delta_Y[indj,indi]%*%TEMP #at each level
   }
  }
  
  #precalc:
  SA = matrix(NA,ncol=r,nrow=I)
  SB = matrix(NA,ncol=r,nrow=I)
  for(l in 1:r) { #for each phi-params.
   deltaGl = deltaG[[l]]
   deltaGlEXP = deltaGl%*%E_delta_Y
   deltaGlVAR = deltaGl%*%Var_delta_Y
   for(i in 1:I) {
    indi <- c(((1:K)-1)*I+i) #indices of all effect at level i
    TEMP = E_delta_Y[indi]%*%SZTZ_i[[i]]
    SUMA1 = 0
    SUMB1 = 0
    for(k in 1:I) {
     indk <- c(((1:K)-1)*I+k) #indices of all effect at level k
     SUMA1 = SUMA1 + sum(SZTZ_i[[i]]*(Var_delta_Y[indi,indk]%*%deltaGlVAR[indk,indi]))
     SUMB1 = SUMB1 + TEMP%*%(Var_delta_Y[indi,indk]%*%deltaGlEXP[indk])
    }
    SA[i,l] = SUMA1
    SB[i,l] = SUMB1
   }
  } 

  COV_dTZTZddTdeltaGld = matrix(0,nrow=I,ncol=r)
  for(l in 1:r) { #for each phi-params.
   deltaGl = deltaG[[l]]
   COV_ddTdeltaGld[,l] = 2*Var_delta_Y%*%deltaGl%*%E_delta_Y
   for(i in 1:I) {
    indi <- c(((1:K)-1)*I+i) #indices of all effect at level i
    COV_dTZTZddTdeltaGld[i,l] = 2*SA[i,l] + 4*SB[i,l]
   }
  }
  
  #BETA
  Var_scorebeta = matrix(0,p,p)
  for(i in 1:I) {
   for(j in 1:I) {
    indi <- c(((1:K)-1)*I+i) #indices of all effect at level i
    indj <- c(((1:K)-1)*I+j) #indices of all effect at level j
    Var_scorebeta = Var_scorebeta + SXTZ_i[[i]]%*%Var_delta_Y[indi,indj]%*%t(SXTZ_i[[j]])/(gamma[i]*gamma[j])
   }
  }
  #make symmetry
  for(i in 1:(p-1))
  for(j in (i+1):p) {
   Var_scorebeta[j,i] = Var_scorebeta[i,j]
  }

  #PHI
  Var_scorephi = matrix(0,r,r)
  for(k in 1:r)
  for(l in 1:r) {
   deltaGi = deltaG[[k]]
   deltaGj = deltaG[[l]]
   deltaGiVar = deltaGi%*%Var_delta_Y
   deltaGjVar = deltaGj%*%Var_delta_Y
   Var_scorephi[k,l] = 1/4*( 2*sum(deltaGiVar*deltaGjVar)+4*t(E_delta_Y)%*%deltaGiVar%*%deltaGj%*%E_delta_Y )
  }
  #Make symmetry
  for(i in 1:(r-1))
  for(j in (i+1):r) {
   Var_scorephi[j,i] = Var_scorephi[i,j]
  }

  #GAMMA
  Var_scoregamma = matrix(0,I,I)
  for(i in 1:I) {
   for(j in i:I) {
    indi <- c(((1:K)-1)*I+i) #indices of all effect at level i
    indj <- c(((1:K)-1)*I+j) #indices of all effect at level j
    COV_didTZTZdj = COV_ddTZTZd[indi,j]
    COV_djdTZTZdi = COV_ddTZTZd[indj,i]
    covval = COV_dTZTZddTZTZd[i,j] + 2*(t(beta)%*%SXTZ_i[[i]]-SYTZ_i[i,])%*%COV_didTZTZdj + 2*(t(beta)%*%SXTZ_i[[j]]-SYTZ_i[j,])%*%COV_djdTZTZdi + 4*(SYTZ_i[i,]-t(beta)%*%SXTZ_i[[i]])%*%Var_delta_Y[indi,indj]%*%t(SYTZ_i[j,]-t(beta)%*%SXTZ_i[[j]])
    Var_scoregamma[i,j] = covval/(4*gamma[i]^2*gamma[j]^2)
   }
  }
  #fill in symmetry:
  for(i in 1:(I-1)) {
   for(j in (i+1):I) {
    Var_scoregamma[j,i] = Var_scoregamma[i,j]
   }
  }

  #GAMMABETA
  Cov_scoregammabeta = matrix(0,nrow=I,ncol=p)
  for(i in 1:I) {
   CovSUM1 = matrix(0,K,p) #Sum_j cov(delta_i,delta_j)tau_jSXTZ_j
   CovSUM2 = matrix(0,1,p)
   for(j in 1:I) {
    indi <- c(((1:K)-1)*I+i) 
    indj <- c(((1:K)-1)*I+j) 
    CovSUM1 = CovSUM1 + Var_delta_Y[indi,indj]%*%t(SXTZ_i[[j]])/gamma[j]
    CovSUM2 = CovSUM2 + t(COV_ddTZTZd[indj,i])%*%t(SXTZ_i[[j]])/gamma[j]
   }
   Cov_scoregammabeta[i,] = (2*(SYTZ_i[i,]-t(beta)%*%SXTZ_i[[i]])%*%CovSUM1 - CovSUM2)/(2*gamma[i]^2)
  }
  #PHIBETA
  Cov_scorephibeta = matrix(0,nrow=r,ncol=p)
  for(l in 1:r) {
   deltaGl = deltaG[[l]]
   CovSUM = rep(0,p)
   for(j in 1:I) {
    indj <- c(((1:K)-1)*I+j) 
    CovSUM = CovSUM + 0.5*t(COV_ddTdeltaGld[indj,l])%*%t(SXTZ_i[[j]])/gamma[j]
   }
   Cov_scorephibeta[l,] = CovSUM
  }
  
  #GAMMAPHI
  Cov_scoregammaphi = matrix(0,nrow=I,ncol=r)
  for(l in 1:r) {
   for(i in 1:I) {
    indi <- c(((1:K)-1)*I+i) #indices of all effect at level i
    Cov_scoregammaphi[i,l] = (2*(SYTZ_i[i,]-t(beta)%*%SXTZ_i[[i]])%*%COV_ddTdeltaGld[indi,l] - COV_dTZTZddTdeltaGld[i,l])/(4*gamma[i]^2)
   }
  } 
  Var_score = rbind(cbind(Var_scorebeta,t(Cov_scorephibeta),t(Cov_scoregammabeta)),
                    cbind(Cov_scorephibeta,Var_scorephi,t(Cov_scoregammaphi)),
                    cbind(Cov_scoregammabeta,Cov_scoregammaphi,Var_scoregamma) )
  hessian = Edeltasql + Var_score
  J = - hessian  #est. observed information matrix

  Var_thetahat = solve(J) 
  stderr = sqrt(diag(Var_thetahat))
  } ) #end stopclock
  show(paste("Timeusage: ",clock(time), sep=""))
  CI = cbind(theta + qnorm(alpha/2)*stderr,theta,theta + qnorm(1-alpha/2)*stderr)
  colnames(CI) = c(paste("(",alpha/2,"%",sep=""),"est",paste(1-alpha/2,"%)",sep=""))
  rownames(CI)=c(paste("beta",1:p,sep=""),paste("phi",1:r,sep=""),paste("gamma",1:length(levelnames),sep=""))
  lmeCI = list(CI=CI,hessian=hessian,Var_thetahat=Var_thetahat,timeusage=time)
  return(lmeCI)
}
