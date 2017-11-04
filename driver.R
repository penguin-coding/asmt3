source('simulation.r')

set.seed(12**3*4**5)

##### WARNING #####
# This simulation takes a long time to run - roughly 1.5 hours on a 4.0 GHz
# processor. The global environment produced by this file is made
# available in the project zip file as 'SimGlobalEnv.r'.

print(Sys.time())
norm.sim <- simulation(dist.func=rnorm,
                       simulations=1000,
                       sample.n=c(10,50,100,500,1000),
                       boot.n=c(9,49,99,199,499,999,1999),
                       boot.method=c('percentile','BCa','parametric','smooth'),
                       stat.func=mean,
                       smooth.sd=0.1)

pois.sim <- simulation(dist.func=rpois,
                       simulations=1000,
                       sample.n=c(10,50,100,500,1000),
                       boot.n=c(9,49,99,199,499,999,1999),
                       boot.method=c('percentile','BCa','parametric','smooth'),
                       stat.func=mean,
                       smooth.sd=0.1,
                       lambda=100)

print(Sys.time())