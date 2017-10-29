#Set of functions students can use in MT4113 Assignemnt 3
#Author: Len Thomas
#Last updated: 19th Oct 2017

#-----------------------------------------------------------------------------------

get.zhat0<-function(est,boot.est){
#Purpose: Return the bias correction factor, zhat0, in BCa bootstrap CI method
#Inputs:
# est - estimated quantity of interest from data
# boot.est - vector of bootstrap estimates of quantity of interest
  
  prop.less<-sum(boot.est<est)/length(boot.est)
  zhat<-qnorm(prop.less)
  return(zhat)
}

#-----------------------------------------------------------------------------------

get.ahat<-function(data,est,fun,...){
#Purpose: Return the acceleration factor, ahat, in BCa bootstrap CI method
#Inputs:
# data - vector of data
# est - estimated quantity of interest from data
# fun - function that can be used to produce est from data via fun(data)
# ... - any other arguments that need to be passed to fun
#Implementation note:
# 1. The routine calls fun(data,...) and expects it to return a scalar equal to est
# 2. Requires a data vector of length at least 2
      
  #Check data vector length
  n<-length(data)
  if(n<2) stop("data vector must be at least length 2\n")

  #Get jacknife estimates of quantity of interest
  jack.est<-numeric(n)
  for(i in 1:n){
    jack.data<-data[-i]
    jack.est[i]<-fun(jack.data)
  }
  
  #Compute ahat
  mean.jack.est<-mean(jack.est)
  ahat.numerator<- sum((mean.jack.est-jack.est)^3)
  ahat.denominator<-6*(sum((mean.jack.est-jack.est)^2))^1.5
  ahat<-ahat.numerator/ahat.denominator
  return(ahat)
}

#-----------------------------------------------------------------------------------
