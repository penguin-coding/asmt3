source('simulation.r')

set.seed(12**3*4**5)

##### WARNING #####
# This simulation takes a long time to run - roughly 1.5 hours on a 4.0 GHz
# processor. The global environment produced by this file is made
# available in the project zip file as 'SimGlobalEnv.r'.

print(Sys.time())

# Full simulation with normal data
norm.sim <- simulation(dist.func=rnorm,
                       simulations=1000,
                       sample.n=c(10,50,100,500,1000),
                       boot.n=c(49,199,499,999,1999,2999),
                       boot.method=c('percentile','BCa','parametric','smooth'),
                       stat.func=mean,
                       smooth.sd=0.1)

# Full simulation with Poisson data
pois.sim <- simulation(dist.func=rpois,
                       simulations=1000,
                       sample.n=c(10,50,100,500,1000),
                       boot.n=c(49,199,499,999,1999,2999),
                       boot.method=c('percentile','BCa','parametric','smooth'),
                       stat.func=mean,
                       smooth.sd=0.1,
                       lambda=100)

# Partial simulation with gamma data, to look at failure tendency with an
# extremely skewed distribution. Exploratory analysis of the above data
# seemed to indicate there was no obvious relationship between sample size and
# failure tendency or between bootstrap resamples and failure tendency.
gamm.sim <- simulation(dist.func=rgamma,
                       simulations=1000,
                       sample.n=200,
                       boot.n=c(99,499,999),
                       boot.method=c('percentile','BCa','smooth'),
                       stat.func=mean,
                       smooth.sd=0.1,
                       shape=3,
                       rate=10)

smooth.values <- seq(0,1,length=30)
smooth.sim.coverage = smooth.sim.length <- rep(NA, length(smooth.values))


print(Sys.time())

norm.results <- calculate.summaries(norm.sim, 0)
pois.results <- calculate.summaries(pois.sim, 100)
gamm.results <- calculate.summaries(gamm.sim, 3/10)

for (results in list(norm.results, pois.results)){
  for (statistic in c('coverage','length','failure tendency')){
    plot(results, statistic=statistic)
  }
}

for (method in 1:4){
  for (statistic in c('failure tendency')){
    
    # Because parametric bootstraps have a coverage of 1 in these simulations,
    # the failure tendency is NaN, since there are no failures. Hence we skip
    # this case to avoid errors with plot.sim.3D:
    if (statistic!='length' & method==3) next
    
    sim.plot.3D(norm.results, statistic=statistic, method=method, hist=F)
    
  }
}

