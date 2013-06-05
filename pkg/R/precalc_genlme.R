#' @title precalc_genlme
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description precalc_genlme calculates sufficient statistics from data which is further used in the genlme-function.
#' @export
#' @details For large amount of data it will be useful to reduce the information into sufficience statistics used for making inference with Maximum Likelihood estimator.
#' 
#' @param dat Dataset with responses, categorize (group), covariates and random effects variables. List elements: Y is Response variable, F is Categorize (name of group, X is Covariates belonging to fixed effects, Z is Covariates belonging to random effects.
#' @param levelnames Name of levels in the one-way categories.
#' @return List of precalculated sufficience matrices needed as input for the genlme function.
#' @references Master Thesis Oyvind Bleka
#' @keywords LME,EM

precalc_genlme <-
function(dat,levelnames) {

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

 #arranged order are as for "levelnames"
 time = system.time( {
 p = dim(dat$X)[2]
 K = dim(dat$Z)[2]
 I = length(levelnames)
 sumsq = function(x){ sum(x^2)} #returns sumsquares
 agg = aggregate(dat$Y,by=list(dat$F),length)
 agglvls = agg[,1]
 gind = compareIndices(levelnames,agglvls) #make agg.ed data in right order.
 #agglvls[gind]==levelnames #OK. Now use index on all agg.ed-data
 n_i = agg[gind,2]
 SYTY_i = aggregate(dat$Y,by=list(dat$F),sumsq)[gind,2]
 SYTX_i = as.matrix(aggregate(dat$Y*dat$X,by=list(dat$F),sum)[gind,2:(1+p)]) #(Ixp) matrix
 SYTZ_i = as.matrix(aggregate(dat$Y*dat$Z,by=list(dat$F),sum)[gind,2:(1+K)]) #(IxK) matrix
 SXTX_i = matrix(list(),ncol=1,nrow=I) #list with (pxp) matrices 
 SZTZ_i = matrix(list(),ncol=1,nrow=I) #list with (pxp) matrices 
 SXTZ_i = matrix(list(),ncol=1,nrow=I) #list with (pxK) matrices 
 for(l in 1:I) {
  S = matrix(0,ncol=p,nrow=p)
  subX = matrix(dat$X[dat$F==levelnames[l],],ncol=p) #take out subset of current level
  subZ = matrix(dat$Z[dat$F==levelnames[l],],ncol=K) #take out subset of current level
  SXTX_i[[l]] = t(subX)%*%subX
  SZTZ_i[[l]] = t(subZ)%*%subZ
  SXTZ_i[[l]] = t(subX)%*%subZ
 }
 } )
 show(paste("Timeusage: ",clock(time), sep=""))
 return(list(SYTY_i=SYTY_i,SXTX_i=SXTX_i, SZTZ_i=SZTZ_i,SYTX_i=SYTX_i,SXTZ_i=SXTZ_i,SYTZ_i=SYTZ_i,levelnames=levelnames,n_i=n_i,timeusage=time))
}
