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

non.par.bootstrap <- function(data, n=999, alpha = 0.05, func = mean,
                           method = 'percentile', forlp = FALSE){
  # purpose : produces a 1 - alpha % confidence interval for the 'func' of the
  #           data, using a non-parametric bootstrap of n samples of size 'size'
  #
  # inputs  : data   - numeric vector of observations from a univariate 
  #                    distribution
  #           n      - the number of resamples to perform for the bootstrap
  #           alpha  - a number between 0 and 1, such that a 1-alpha %
  #                    confidence interval will be produced
  #           func   - function to be used to calculate the statistic of
  #                    interest for each sample. The default is 'mean'. 
  #           method - 'percentile' or 'bca'. Specifies the method to be used to 
  #                    obtain the quantiles to sample from to obtain the
  #                    interval
  #           forlp  - default is FALSE. When TRUE, the function uses a for 
  #                    loop to conduct the bootstrap. This functionality was 
  #                    added in so that we could verify that the non-loop 
  #                    version was in fact faster. 
  #          
  # output  : named vector containing the lower and upper bounds of the interval
  
  # Input checks and setting default values: 
  if (class(data)!='numeric' & class(data)!='integer'){
    stop('input data must be numeric')}

  if (n%%1!=0 | n<1) stop('n must be a positive integer')
  if (alpha<0 | alpha>1) stop('alpha must be between 0 and 1')
  func <- match.fun(func) # to allow the user to pass in the func or func name
  if (!is.function(func)) stop('invalid function supplied as func argument')
  if (!method %in% c('percentile','BCa')) stop('invalid method')
  
  samples <- non.parametric.sample(data, n)       # generate the random samples
  samples <- cbind(samples,data)                  # add in the observed data
  samples <- apply(samples,2,func)                # calculate statistics
  
  lower <- alpha/2      # percentile method
  upper <- 1 - alpha/2  # intervals
  
  if (method=='BCa'){
    alphas <- get.bca.alphas(data, func(data), samples, alpha, func)
    lower <- alphas[1]
    upper <- alphas[2]
    }
  
  CI <- quantile(samples, probs=c(lower, upper))
  return(CI)
}

simulation <-  function(dist.func, simulations, sample.n, boot.n, boot.method,
                        stat.func, alpha, ...){
  # purpose : run a set of simulations
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
  #           boot.method - 'percentile' or 'BCa' as a character input. Can be
  #                         a vector
  #           stat.func   - The function which calculates the statistic of
  #                         interest for which we are producing a confidence 
  #                         interval for
  #           alpha       - We produce (1-alpha)*100 % confidence intervals
  #           ...         - Extra parameters to be passed to dist.func
  # 
  # output : a multi-dimensional array with named dimensions, containing all of 
  #          the statistics produced by the simulated intervals. Has class 
  #          'simulation.output.object'
  
  # generate the multi-dimensional array which will store all of the generated 
  output <- array(data = NA, # intervals
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
          boot <-  non.par.bootstrap(dataset,n=boot.n.setting,
                                    alpha=alpha,func = stat.func,
                                    method = boot.method.setting)
          
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
