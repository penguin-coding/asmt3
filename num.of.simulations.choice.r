source('simulation.r')
dosims <- function(x,sample.n,boot.n){
    
    norm.sim <- simulation(dist.func=rnorm,
                           simulations=x,
                           sample.n=sample.n,
                           boot.n=boot.n,
                           boot.method=
                             c('percentile','BCa','parametric','smooth'),
                           stat.func=mean,
                           smooth.sd=0.1)
  
    pois.sim <- simulation(dist.func=rpois,
                           simulations=x,
                           sample.n=sample.n,
                           boot.n=boot.n,
                           boot.method=
                             c('percentile','BCa','parametric','smooth'),
                           stat.func=mean,
                           smooth.sd=0.1,
                           lambda=100)
  
    gamm.sim <- simulation(dist.func=rgamma,
                           simulations=x,
                           sample.n=sample.n,
                           boot.n=boot.n,
                           boot.method=c('percentile','BCa','smooth'),
                           stat.func=mean,
                           smooth.sd=0.1,
                           shape=3,
                           rate=10)
  
  smooth.sd <- seq(0.01,0.25,length=10)
  smoot.sim <- as.list(rep(NA, length(smooth.sd)))
  
  counter = 0
  for (i in smooth.sd){
    counter = counter + 1
    sim <- simulation(dist.func=rnorm,
                      simulations=x,
                      sample.n=500,
                      boot.n=999,
                      boot.method='smooth',
                      stat.func=mean,
                      smooth.sd=smooth.sd[counter])
    smoot.sim[[counter]] <- sim[1,1,1,]
  }
  
  
  return()
}

get.time = function(x,sample.n,boot.n){
  return(system.time(dosims(x,sample.n,boot.n))[3])
}

x <- seq(2)
y <-  lapply(x, get.time, sample.n=c(50,100,500,1000), boot.n=c(99,199,499,999))
y <- unlist(y)

A = lm(y~x-1)  # We fit a lm with no intercept, and check
summary(A)     # it produces a reasonable answer. 

plot(x,y)                    # We visualise the fitted line
abline(0,A$coefficients)

pred.time <- function(s){     # We write a very simple function
  secs <- A$coefficients*s    # which returns the expected
  mins <- secs/60             # calculation time for s simulations
  hours <- mins/60            # per case, based on our observed
  return(hours)               # values
}

plot.x = seq(1,2000)        # we plot expected calculation times
plot.y = pred.time(plot.x)  # so that we can choose how many to use

plot(plot.x, plot.y, xlab = 'number of simulations per case',
     ylab = 'Time for all simulations to run - hours',
     type='l',
     main = 'Predicted calculation times for 1.51GHz CPU')