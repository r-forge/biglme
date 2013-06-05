#' @title genlmesim
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description genlmesim simulates dataset from a Linear Mixed Effects model which may assume any general covariance structures for the effects.
#' @export
#' @details The model assumes gaussian noise with different variance paramater for each group (heterogenous variance structure) and gaussian random effect with any covariance strucuture.
#' @param model A list with function-elements' G' and' ldetG' taking a' phi'-parameter.
#' @param true A list with vector-elements' beta','phi' and' gamma' used as parameter in the LME model.
#' @param designM A list with covariation matrix' X' belonging to fixed effects and random effect matrix' Z' belonging to random effects.
#' @param n_i Number of data for each categorizations.
#' @param levelname Name of levels in the one-way categories.
#' @return Dataset with responses, categorize (group), covariates and random effects variables
#' \item{Y}{Response variable.}
#' \item{F}{Categorize (name of group).}
#' \item{X}{Covariates belonging to fixed effects.}
#' \item{Z}{Covariates belonging to random effects.}
#' \item{delta}{Sampled random effects.}
#' @keywords LME,EM

genlmesim <-
function(model=NULL,true=NULL,designM=NULL, n_i=NULL,levelname=NULL) {
 if(is.null(model) | is.null(true) | is.null(n_i)) stop("Missing arguments")
 G <- model$G(true$phi)
 q = dim(G)[1]
 I = length(n_i) #number of levels
 Z = as.matrix(designM$Z)
 K = dim(Z)[2]
 X = as.matrix(designM$X)
 beta <- true$beta
 gamma <- true$gamma
 n = sum(n_i)
 
 #simulate delta's first
 delta = chol(solve(G))%*%rnorm(q)
 y = rep(NA,n)
 F = rep(NA,n)
 cn_i = c(0,cumsum(n_i))+1 #start indices for each levels
 for(i in 1:I) {
  ind = c(((1:K)-1)*I+i) #indices having same levels
  range = cn_i[i]:(cn_i[(i+1)]-1) #indice range of level i
  y[range] = c(Z[range,]%*%as.matrix(delta[ind])) + c(X[range,]%*%as.matrix(beta))  + sqrt(gamma[i])*rnorm(n_i[i])
  F[range] = i
  if(!is.null(levelname)) F[range] = levelname[i]
 }
 return(list(Y=y,F=F,X=X,Z=Z,delta=delta))
}
