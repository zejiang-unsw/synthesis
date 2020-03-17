library(MASS)
#----------------------------------------------------------------------------------------------------
data.gen.nl<-function(nobs,ndim=9,noise=1){
  #nobs<-1000;ndim=9;noise=0
  nwarm=500
  n=nobs+nwarm
  x<-matrix(NA,n,1)
  dp<-matrix(NA,n,ndim) 
  
  for(i in 1:ndim) dp[,i]<-rnorm(n,mean=0,sd=1)
  #plot.ts(dp[,1:10])
  
  for (i in 1:n){
    eps<-rnorm(1,mean=0,sd=1)
    x[i]<- 0.8*dp[i,1]^2 + 0.4*dp[i,2] + 0.7*dp[i,3] + noise*eps
  }
  
  x=x[(nwarm+1):n]
  dp=dp[(nwarm+1):n,]
  data_generated<-list(x=x,dp=dp,true.cpy=c(1,2,3))
  #plot.ts(cbind(x,dp[,c(2,6,9)]),type="l")
  
  return(data_generated)
}	
#----------------------------------------------------------------------------------------------------
data.gen.nl1<-function(nobs,ndim=9,noise=0){
  #nobs<-1000;ndim=9;noise=0
  nwarm=500
  n=nobs+nwarm
  x<-matrix(NA,n,1)
  dp<-matrix(NA,n,ndim) 
  
  for(i in 1:ndim) dp[,i]<-rnorm(n,mean=0,sd=1)
  #plot.ts(dp[,1:10])
  
  for (i in 1:n){
    eps<-rnorm(1,mean=0,sd=1)
    x[i]<- dp[i,2]^2 + cos(dp[i,6]) + 0.35*sin(dp[i,9]) + noise*eps
  }
  
  x=x[(nwarm+1):n]
  dp=dp[(nwarm+1):n,]
  data_generated<-list(x=x,dp=dp,true.cpy=c(2,6,9))
  #plot.ts(cbind(x,dp[,c(2,6,9)]),type="l")

  return(data_generated)

}	
#----------------------------------------------------------------------------------------------------
data.gen.nl2<-function(nobs,ndim=9,r=0.6,noise=1){
  #nobs<-1000;ndim=9;noise=1;r=0.6
  nwarm=500
  n=nobs+nwarm
  x<-matrix(NA,n,1)
  dp<-matrix(NA,n,ndim) 
  
  Sigma=matrix(rep(r,ndim*ndim),ndim); diag(Sigma) <- 1
  dp <-mvrnorm(n,mu=rep(0,ndim),Sigma=Sigma,  empirical=TRUE)
  #plot.ts(dp[,1:10])
  
  for (i in 1:n){
    eps<-rnorm(1,mean=0,sd=1)
    x[i]<- dp[i,2]^2 + cos(dp[i,6]) + 0.35*sin(dp[i,9]) + noise*eps
  }
  
  x=x[(nwarm+1):n]
  dp=dp[(nwarm+1):n,]
  data_generated<-list(x=x,dp=dp,true.cpy=c(2,6,9))
  #plot.ts(cbind(x,dp[,c(2,6,9)]),type="l")
  
  return(data_generated)
}	



