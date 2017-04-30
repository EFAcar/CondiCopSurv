#' @title Simulating Censoring Times for GLRT Bootstrap Samples
#' @description This function simulates censoring times in boostrap samples of right censored survival data
#'
#' @param Y1,Y2 vectors of event times
#' @param status1,status2 vectors of censoring indicators
#' @importFrom survival survfit
#' @keywords internal
BootSim.Censor = function(Y1,Y2, status1, status2){
  
  n = length(Y1)
  nevent_cluster = status1 + status2
  max_cluster = pmax(Y1,Y2)
  statusmax_cluster=status1 * status2
  if(sum(statusmax_cluster) == n){cstar = rep(Inf,n)
  }else{
    
    SKMcens=summary(survival::survfit(Surv(max_cluster,1-statusmax_cluster)~1))
    ncens=length(SKMcens$time)
    
    probcond=matrix(0,nrow=n , ncol=ncens) ; cstar=matrix(NA,nrow=n)
    
    for(i in 1:n){
      if(nevent_cluster[i]==2 && max_cluster[i] < SKMcens$time[ncens]){
        if(max_cluster[i]<SKMcens$time[1]){
          probcond[i,1]=1-SKMcens$surv[1]
          for(j in 2:ncens){
            probcond[i,j]=SKMcens$surv[j-1]-SKMcens$surv[j]
          }}
        
        for(w in 2:ncens){
          if(SKMcens$time[w-1]<=max_cluster[i] && max_cluster[i]<SKMcens$time[w]){
            for(j in w:ncens){
              probcond[i,j]=(SKMcens$surv[j-1]-SKMcens$surv[j])/SKMcens$surv[w-1]
            }}}
        cstar[i]=sample(SKMcens$time,1,prob=probcond[i,])
      }
      if(nevent_cluster[i]==2 && max_cluster[i]>=SKMcens$time[ncens]){cstar[i]=max_cluster[i]}
      if(nevent_cluster[i]!=2){cstar[i]=max_cluster[i]}
    }}
  
  return(cstar)
}
