# I confirm that the attached is my own work, except where clearly indicated in
# the text.

source('simulation.r')
library('magrittr')   # ceci n'est pas une %>%

set.seed(12**3*4**5)  # to make results reproducible to the reader

##### WARNING #####

# This simulation takes a long time to run - roughly 3 hours on a 4.0 GHz
# processor. The global environment produced by this file is made
# available in the project zip file as 'SimGlobalEnv.r'.

print(Sys.time()) # We also print the end time

sample.n <- c(50,100,500,1000)
boot.n <- c(99,199,499,999)
simulations <- 1000
boot.method=c('percentile','BCa','parametric','smooth','par.fit')
# Full simulation with normal data
norm.sim <- simulation(dist.func='rnorm',
                       simulations=simulations,
                       sample.n=sample.n,
                       boot.n=boot.n,
                       boot.method=boot.method,
                       stat.func=mean,
                       smooth.sd=0.1,
                       mean=0,
                       sd= 1)

# Full simulation with Poisson data
pois.sim <- simulation(dist.func='rpois',
                       simulations=simulations,
                       sample.n=sample.n,
                       boot.n=boot.n,
                       boot.method=boot.method,
                       stat.func=mean,
                       smooth.sd=0.1,
                       lambda=100)

# Full simulation with gamma data. The skewness of this distribution will 
# hopefully allow for more meaningful insight into the behaviour of failure
# tendency for each type of bootstrap. We exempt parametric bootstraps, 
# because their coverage for high numbers of bootstrap resamples is 1.
gamm.sim <- simulation(dist.func='rgamma',
                       simulations=simulations,
                       sample.n=sample.n,
                       boot.n=boot.n,
                       boot.method=boot.method,
                       stat.func=mean,
                       smooth.sd=0.1,
                       shape=3,
                       rate=10)





# Side analysis looking at how the quantity of noise affects smooth bootstraps #

# We can't vary the proportion of the variance of the smoothing terms
# with respect to the observed data using the simulation function, so 
# we manually produce a small simulation dataset for analysis:
smooth.sd <- seq(0.01,0.33,length=20)            # set our chosen values
smoot.sim <- as.list(rep(NA, length(smooth.sd))) # object to store outputs

counter = 0
for (i in smooth.sd){
  counter = counter + 1
  sim <- simulation(dist.func=rnorm,
                    simulations=simulations,
                    sample.n=500, # we fix sample.n and boot.n for this
                    boot.n=999,   # simulation, to save on computing time.
                    boot.method='smooth',
                    stat.func=mean,
                    smooth.sd=smooth.sd[counter])
  smoot.sim[[counter]] <- sim[1,1,1,]
}
remove(i) ; remove(counter) # remove unnecessary global vars
print(Sys.time())





### side analysis looking at BCa vs percentile for low n and high resamples ####

pois.sim2 <- simulation(dist.func='rpois',
                        simulations=simulations,
                        sample.n=c(10,20,35,50),
                        boot.n=c(2999,3499),
                        boot.method=c('percentile','BCa'),
                        stat.func=mean,
                        lambda=100)





############# Producing plots to visualise the simulation results ##############

norm.results <- calculate.summaries(norm.sim, 0)    # calculate coverage, length
pois.results <- calculate.summaries(pois.sim, 100)  # and 'failure tendency'
gamm.results <- calculate.summaries(gamm.sim, 3/10) # for each simulation

for (results in list(norm.results, pois.results, gamm.results)){
  
  for (statistic in c('coverage','length','failure tendency')){
    
    # to avoid "the condition has length > 1" warnings, we compare the first
    # element only - which should be sufficient for our purposes:
    if (results[1,1,1,1]==norm.results[1,1,1,1]){main <- 'normal'}
    else if (results[1,1,1,1]==pois.results[1,1,1,1]){main <- 'poisson'}
    else {main <- 'gamma'}
    
    main <- paste(main, 'simulated deviates')
    
    plot(results, statistic=statistic, main=main)
  }
}

smooth.stats <- array(dim=c(length(smoot.sim),3))

for (i in 1:length(smoot.sim)){                         # we manually create a
  cov <- get.coverage(smoot.sim[[i]],0)                 # a matrix of summary 
  len <- get.length(smoot.sim[[i]])                     # statistics for the
  f.t <- get.coverage(smoot.sim[[i]],0,failure.t = T)   # simulation where only
  smooth.stats[i,] <- c(cov,len,f.t)                    # the smooth.sd
}                                                       # quantity varies


# We manually produce some plots for the smooth bootstrap analysis:
for (i in 1:3){
  ylab <-  switch(i, '1'='coverage', '2'='length', '3'='failure tendency')
  main <- 'Smooth bootstrap'
  y <- smooth.stats[,i]
  plot(smooth.sd, y,xlab = 'smooth.sd', ylab=ylab, main=main,
       col=4)
}


# Produce the plot which supports our claim that the BCa is better for low n:
pois.sim2 %>% calculate.summaries(100) %>% plot(statistic='coverage',
                                                main='Poisson data')





################ Create some 3D plots for the standard analyses ################

for (method in 1:4){
  for (statistic in c('coverage')){

    # Because parametric bootstraps have a coverage of 1 in these simulations,
    # the failure tendency is NaN, since there are no failures. Hence we skip
    # this case to avoid errors with plot.sim.3D:
    if (statistic!='length' & method==3) next

    sim.plot.3D(norm.results, statistic=statistic, method=method, hist=F)

  }
}