library(boot)
source('simulation.r')

# given the strange BCa results, I got very worried that my BCa code was
# wrong. Not finding an issue with the formulas being used to retrieve the
# alpha values, I decided a comparison with the boot package's own BCa would be
# the next reasonable point of call as a means of testing. The following is a
# small piece of code which checks that my BCa and the boot package BCa perform
# similarly in a small handful of situations (for poisson data).

our.mean <-  function(d,w){return(mean(d[w]))} # needed by boot function

compare.BCAs <- function(sims, n, R, statfuncs){
  # purpose : Compare performance of the boot package BCa bootstrap and 
  #           The one I wrote for this project
  # 
  # inputs  : sims      - the number of simulations to run for the test case
  #           n         - sample size for observations
  #           R         - number of bootstrap resamples
  #           statfuncs - list of functions which calculate relevant 
  #                       statistics. The first entry should be a function 
  #                       of the type required for my BCa, the second should 
  #                       be of the type required by the boot package's
  #                       boot function
  #
  # output  : The observed coverage of each method, as a percentage.
  
  MyBCa <- vector(length=sims*2)   # Create blank vectors to store
  BootBCa <- vector(length=sims*2) # the bootstrap intervals
  
  for (i in seq(1,2*sims,2)){
    
    data = rpois(n, lambda=100)  # generate the data
    
    A <-  bootstrap(data, n=R, func=statfuncs[[1]],           # generate my own
                               method='BCa',check.inputs=F)   # BCa bootstrap
    
    resamps <- boot(data, statfuncs[[2]], R)                  # generate the
    interval <- boot.ci(resamps, type='bca')$bca[4:5]         # boot interval
    
    MyBCa[i] <- A[1]              # add my interval
    MyBCa[i+1] <- A[2]            # to the output
    
    BootBCa[i] <- interval[1]     # add the boot interval
    BootBCa[i+1] <- interval[2]   # to the output
  }
  
  coverages <- c(get.coverage(MyBCa, 100), get.coverage(BootBCa,100))
  names(coverages) <- c('mine', 'boot.package')
  return(coverages)
}

set.seed(666)
R <- c(999)
N <- c(20)
sims <- 10000

system.time(
for (r in R){
  for (n in N){
    cat('\n Coverages for sample size',n,'and',r,'resamples:\n')
    print(compare.BCAs(sims, n, r, list(mean, our.mean)))
  }
}
)[3]


# Playing around with different settings in the above lines of code reveals
# that my method does in fact have a consistent small drop in percentage
# cover, but I still have no idea why.