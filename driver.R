source('simulation.r')

set.seed(18061996)

norm.sim <- simulation(dist.func=rnorm,
                       simulations=500,
                       sample.n=c(10,50,100,500,1000),
                       boot.n=c(9,19,29,99,199,499,999),
                       boot.method=c('percentile','BCa'),
                       stat.func=mean,
                       alpha=0.05)

poisson.sim <- simulation(dist.func=rpois,
                          simulations=500,
                          sample.n=c(10,50,100,500,1000),
                          boot.n=c(9,19,29,99,199,499,999),
                          boot.method=c('percentile','BCa'),
                          stat.func=mean,
                          alpha=0.05,
                          lambda=10)

smal.nor <- simulation(dist.func=rnorm,
                       simulations=100,
                       sample.n=c(10,50,100),
                       boot.n=c(49,99,199),
                       boot.method=c('percentile','BCa'),
                       stat.func=mean,
                       alpha=0.05)

smal.poi <- simulation(dist.func=rpois,
                       simulations=100,
                       sample.n=c(10,50,100),
                       boot.n=c(49,99,199),
                       boot.method=c('percentile','BCa'),
                       stat.func=mean,
                       alpha=0.05,
                       lambda=10)
