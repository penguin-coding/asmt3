# I confirm that the attached is my own work, except where clearly indicated
# in the text. 

source('BCaHelperFunctions.r') # Use Len's code to help with BCa Bootstrap

non.parametric.sample <- function(data, n){
  # purpose : produces n random samples of size 'size' from the supplied data
  #
  # inputs  : data - numeric vector of univariate observations
  #           n    - positive integer number of samples to be drawn
  #
  # output  : size*n dimension matrix. Each column is a generated sample.
  
  if (class(data)!='numeric' & class(data)!='integer'){
    stop('input data must be numeric')}
  
  if (n%%1!=0 | n<1) stop('n must be a positive integer')
  
  return(replicate(n, sample(x=data, size=length(data),replace=T)))
}

get.bca.alphas <- function(data, est, boot.est, alpha, func){
  # purpose : calculates values alpha1 and alpha2, to be taken as the quantiles
  #           used to compute a BCa confidence interval
  #
  # inputs  : data     - numeric vector of univariate observations 
  #           est      - estimate of statistic of interest from original sample
  #           boot.est - estimate of statistic of interest from each of the 
  #                      bootstrap samples
  #           alpha    - a value between 0 and 1 such that a (1-a)*100 %
  #                      confidence interval will be produced.
  #           func     - a function such that func(data) produces est. It is the
  #                      function which provides an estimate of the statistic of
  #                      interest given an input dataset. 
  #
  # output  : a named vector of alphas such that quantile(data, probs=alphas)
  #           returns the BCa confidence interval for the statistic of interest
  #
  # notes: - this method is treated as private, and so does not type check its
  #          inputs. It should only be called by functions which have already
  #          checked their inputs and produce consistent outputs of the correct
  #          format.
  
  alpha <-  alpha/2 # transform alpha s.t. we get a (1-2a)*100 % CI
  zhat0 <- get.zhat0(est, boot.est)
  ahat <- get.ahat(data, est, func)
  
  # Calculate the values of alpha1 and alpha2 according to the formula specified
  # in the assignment 3 outline pdf:
  alpha1 <- zhat0 + ( (zhat0 + qnorm(alpha)) / (1-ahat*(zhat0 + qnorm(alpha))) )
  alpha2 <- zhat0 + ((zhat0 + qnorm(1-alpha))/(1-ahat*(zhat0 + qnorm(1-alpha))))
  
  alpha1 <- pnorm(alpha1) ; alpha2 <- pnorm(alpha2)
  alphas <- c(alpha1,alpha2) ; names(alphas) <- c('alpha1','alpha2')
  
  return(alphas)
}

bootstrap <- function(data, n=999, alpha = 0.05, func = mean,
                      method = 'percentile', smooth.sd = 0.2, 
                      dist.func = NULL, ...){
  # purpose : produces a 1 - alpha % confidence interval for the 'func' of the
  #           data, using a non-parametric bootstrap of n samples of size 'size'
  #
  # inputs  : data       - numeric vector of observations from a univariate 
  #                        distribution
  #           n          - the number of resamples to perform for the bootstrap
  #           alpha      - a number between 0 and 1, such that a 1-alpha %
  #                        confidence interval will be produced
  #           func       - function to be used to calculate the statistic of
  #                        interest for each sample. The default is 'mean'. 
  #           method     - 'percentile', 'BCa', 'parametric' or 'smooth'.
  #                        Specifies the bootstrap method to be used. When
  #                        'parametric' is chosen, the percentile method is used
  #                        to calculate the confidence interval using the
  #                        bootstrap samples, and a function from which to
  #                        sample the data must be specified. All remaining 
  #                        options produce non-parametric bootstraps. Option
  #                        'smooth' adds a normal noise centred at 0. The
  #                        chosen standard deviation is a fraction of the sample
  #                        sd. This is set using the parameter 'smooth.sd'.
  #           smooth.sd  - Multiplier for the sample standard deviation. When 
  #                        method = 'smooth', a normal noise id added to 
  #                        each bootstrap resample. It has mean 0 and standard
  #                        deviation smooth.sd * sd(data).
  #           dist.func  - function to sample the data from when parametric is 
  #                        set to TRUE. It is assumed that the first argument
  #                        in any call to dist.func is the number of random
  #                        deviates to be produced, as is the convention with
  #                        rnorm, runif, rpois, rgamma, etc.
  #           ...        - extra optional parameters to be passed to dist.func
  #          
  # output  : named vector containing the lower and upper bounds of the interval
  
  # Input checks and setting default values: 
  if(!is.null(dist.func)) dist.func <- match.fun(dist.func)
  if (class(data)!='numeric' & class(data)!='integer'){
    stop('input data must be numeric')}
  
  if ((method=='parametric' & is.null(dist.func))|
      (method=='parametric' & !is.function(dist.func))){
    stop('when parametric is set to TRUE a dist.func function must be provided')
  }
  
  if (method=='smooth' & (class(smooth.sd)!='numeric' | smooth.sd<0 |
                          length(smooth.sd)>1 )){
    
    stop('When method = \'smooth\', smooth.sd must be a positive scalar')
  }
  
  if (n%%1!=0 | n<2) stop('n must be a positive integer greater than 1')
  if (alpha<0 | alpha>1) stop('alpha must be between 0 and 1')
  func <- match.fun(func) # to allow the user to pass in the func or func name
  if (!is.function(func)) stop('invalid function supplied as func argument')
  if (!method %in% c('percentile','BCa','parametric','smooth')){
    stop('invalid method')}
  
  ### End of input-checks
  if (smooth.sd==0 & method=='smooth') method <- 'percentile'
  
  if (method!='parametric'){  # generate the random samples with replacement
    samples <- replicate(n, sample(x=data, size=length(data),replace=T)) 
    
    if (method=='smooth'){ # add noise to the data for a smooth bootstrap
      noise <- replicate(n, rnorm(length(data), sd = sd(data)*smooth.sd))
      samples <- samples + noise
    }
    
  }
  
  else{ # parametric resamples
    # the replicate function worked badly with functions like rpois and rt, 
    # so a less efficient method has to be used to produce parametric samples:
    samples <- matrix(nrow=length(data),ncol=n)
    for (i in 1:n){samples[,i] <- dist.func(length(data),...)}
  }
  
  samples <- cbind(samples,data)                # add in the observed data
  samples <- apply(samples,2,func)              # calculate statistics
  lower <- alpha/2      # percentile method
  upper <- 1 - alpha/2  # intervals
  
  if (method=='BCa'){
    alphas <- get.bca.alphas(data, func(data), samples, alpha, func)
    lower  <- alphas[1]
    upper  <- alphas[2]
    }
  
  CI <- quantile(samples, probs=c(lower, upper))
  return(CI)
}

simulation <-  function(dist.func, simulations, sample.n, boot.n, boot.method,
                        stat.func=mean, alpha=0.05, smooth.sd=0.2,...){
  # purpose : run a set of simulations with different settings.
  #
  # inputs  : dist.func   - The function which should be used to generate the 
  #                         random data used at the start of each simulation
  #           simulations - The number of simulations to run ; how many
  #                         bootstraps should be produced for each setting?
  #           sample.n    - The sample size to be used every time a sample is
  #                         generated using dist.func. Can be a vector.
  #           boot.n      - The number of resamples each bootstrap should
  #                         perform in order to produce its interval. Can be a
  #                         vector
  #           boot.method - 'percentile', 'BCa', 'smooth' or 'parametric' as a
  #                         character input. Can be a vector
  #           stat.func   - The function which calculates the statistic of
  #                         interest for which we are producing a confidence 
  #                         interval for
  #           alpha       - We produce (1-alpha)*100 % confidence intervals
  #           smooth.sd   - What fraction of the sample sd should the sd of the
  #                         noise added to the data have for a smooth bootstrap?
  #           ...         - Extra parameters to be passed to dist.func
  # 
  # output : a multi-dimensional array with named dimensions, containing all of 
  #          the statistics produced by the simulated intervals. Has class 
  #          'simulation.output.object'
  
  ### note:
  ### we don't type check inputs since the bootstrap function will do it for us.
  
  # generate the multi-dimensional array which will store all of the generated 
  output <- array(data = NA,                                       # intervals
                  dim = c(length(sample.n), length(boot.n),length(boot.method),
                          2*simulations),
                  dimnames = list(paste('sample.n:',as.character(sample.n)),
                                  paste('boot.n:',as.character(boot.n)),
                                  paste('boot.method:',boot.method))
                  )
  
  for (sample.n.setting in sample.n){
    for (boot.n.setting in boot.n){
      for (boot.method.setting in boot.method){
        
        sample.n.index <-  which(sample.n==sample.n.setting) # extract indices
        boot.n.index <- which(boot.n==boot.n.setting)        # of settings
        boot.method.index <- which(boot.method==boot.method.setting)
        
        sims <- matrix(nrow=2, ncol=simulations)
        
        for (i in 1:simulations){
          dataset <- dist.func(sample.n.setting, ...) # get the original sample
          
          # get the bootstrap interval for that dataset:
          boot <-  bootstrap(dataset, n=boot.n.setting, alpha=alpha,
                             func = stat.func,
                             method = boot.method.setting,
                             smooth.sd = smooth.sd,
                             dist.func = dist.func, ...)
          
          # add the bootstrap to the matrix of results:
          sims[,i] <- boot
        }

        # add the set of simulated bootstrap intervals to the output array:
        output[sample.n.index, boot.n.index, boot.method.index,] <-  sims
      }
    }
  }
  class(output) <- 'simulation.output.object'
  return(output)
}

calculate.summaries <- function(simulation.output.object, true.value){
  # purpose : takes as input a simulation.output.object and the true.value of 
  #           the statistic for the distribution used to produce the deviates
  #           and calculates some summaries (coverage, length etc.) using the
  #           simulation results contained in the simulation.output.object
  #
  # inputs  : simulation.output.object - the result of calling the function 
  #           'simulation' which is a multi-dimensional array containing 
  #           all the simulated bootstrap intervals at each level of the 
  #           simulation settings
  #
  # output  : a simulation.summaries object, it is simply a
  #           simulation.output.object with the 4th dimension of the array 
  #           representing the various summaries we have calculated for those
  #           simuation settings
  
  if (class(simulation.output.object)!='simulation.output.object'){
    stop('input must be a valid simulation.output.object')
  }
  
  if (class(true.value)!='numeric') stop('true.value must be a real number')
  
  dims <- dim(simulation.output.object)
  
  output <- array(dim = c(dims[1],dims[2],dims[3],2))
  dimnames(output) <- dimnames(simulation.output.object)
  for (i in 1:dims[1]){     # With nested for loops, go through the 
    for(j in 1:dims[2]){    # simulated bootstrap intervals and calculate
      for(k in 1:dims[3]){  # the summary statistics of interest:
        
        boot.ints <- simulation.output.object[i,j,k,]  # extract intervals
        
        coverage <- get.coverage(boot.ints,true.value) # calculate the
        length <- get.length(boot.ints)                # statistics
         
        summaries <- c(coverage,length)                # add them to the 
        names(summaries) <- c('coverage','length')     # output object with
        output[i,j,k,] <- summaries                    # appropriate names
      }
    }
  }
  
  class(output) <- 'simulation.summary.object'
  return(output)
}

get.coverage <- function(bootstrap.results, true.value){
  # purpose : returns the observed coverage, given a vector which contains
  #           a sequence of confidence intervals
  #
  # input   : bootstrap.results - a vector containing bootstrap intervals in the
  #           format c(lower1, upper1, lower2, upper2, etc.)
  #
  # output  : numeric scalar ; the observed coverage given the vector of
  #           bootstrap intervals
  
  if ( class(bootstrap.results)!='numeric' | class(true.value)!='numeric'){
    stop('invalid input')}
  
  n = length(bootstrap.results)
  if (n%%2!=0) stop('input of odd length is not allowed')
  
  lowers <- bootstrap.results[seq(1,n,2)] # we split our intervals into
  uppers <- bootstrap.results[seq(2,n,2)] # vectors of lower and upper bounds
  
  # is the true.value contained in each of our confidence intervals? :
  in.interval <- (true.value>=lowers & true.value<=uppers)
  
  return(sum(in.interval)/(n/2)) # return the observed coverage
}

get.length <- function(bootstrap.results){
  # purpose : returns the observed average interval length, given a vector which
  #           contains a sequence of confidence intervals
  #
  # input   : bootstrap.results - a vector containing bootstrap intervals in the
  #           format c(lower1, upper1, lower2, upper2, etc.)
  #
  # output  : numeric scalar ; the observed average interval length given the
  #           vector of bootstrap intervals
  
  if ( class(bootstrap.results)!='numeric') stop('invalid input')
  
  n = length(bootstrap.results)
  if (n%%2!=0) stop('input of odd length is not allowed')
  
  lowers <- bootstrap.results[seq(1,n,2)] # we split our intervals into
  uppers <- bootstrap.results[seq(2,n,2)] # vectors of lower and upper bounds
  
  return(mean(abs(uppers-lowers)))# return the estimated average interval length
}

plot.simulation.summary.object <- function(simulation.summary.object,
                                           statistic='coverage'){
  # purpose : plots the statistic of interest for a set of simulation
  #           bootstrap confidence intervals, for all levels of 'factor'. Fixes
  #           the other setting values at their highest setting i.e. uses the 
  #           the largest sample size and bootstrap resamples available.
  #
  # inputs  : simulation.summary.object - array of summary statistics for 
  #                                       simulation intervals.
  #           statistic                 - summary statistic of interest, 
  #                                       'coverage', 'length'
  # output  : None, produces a plot. 
  
  if (class(simulation.summary.object)!='simulation.summary.object'){
    stop('invalid input type')}
  
  if ( !(statistic %in% c('coverage','length')) ){
    stop('invalid choice of statistic')}
  
  # fetch summary statistic index:
  ifelse(statistic=='coverage', stat.ind <- 1, stat.ind <- 2) 
  
  dims <- dim(simulation.summary.object) # get maximum index for each level
  
  ### first plot : sample size
  y <- list()
  
  # Extract the sample sizes from the dimension names using some functions
  # which parse through strings: 
  x <- as.numeric(
    gsub('[^0-9]','',dimnames(simulation.summary.object)[[1]]))
  
  # Draw the first plot:
  method.names <- gsub('boot.method: ','',
                       dimnames(simulation.summary.object)[[3]])

  matplot(x,
          simulation.summary.object[,dims[2],,stat.ind], 
          ylab=statistic, 
          xlab='sample size',
          type='l',
          col=seq(1,dims[3]))
  
  legend('topright',
         method.names, 
         lty=1, col=seq(1,dims[3]), bty='n', cex=.75)
  
  
}

