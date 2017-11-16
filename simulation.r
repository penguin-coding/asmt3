# I confirm that the attached is my own work, except where clearly indicated
# in the text. 

source('BCaHelperFunctions.r') # Use Len's code to help with BCa Bootstrap
library(plot3D)                # Used by the sim.plot.3D function
library(reshape2)              # Used by the sim.plot.3D function
library(rgl)                   # Used by the sim.plot.3D function
library(magrittr) # ceci n'est pas une pipe, hon hon hon

non.parametric.sample <- function(data, n){
  # purpose : produces n random samples of size length(data) from the supplied
  #           data
  #
  # inputs  : data - numeric vector of univariate observations
  #           n    - positive integer number of samples to be drawn
  #
  # output  : length(data)*n dimension matrix. Each column is a generated sample
  
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

bootstrap.type.checks <- function(data, n, alpha, func, method,
                                  smooth.sd, dist.func=NULL){
  # purpose : checks all the inputs of the bootstrap function are of the
  #           expected type and satisfy all required conditions they impose on
  #           each other
  #
  # inputs  : data      - vector of univariate observations
  #           n         - number of bootstrap resamples
  #           alpha     - target coverage of interval
  #           func      - function which calculates statistic of interest
  #           method    - character name of method to be used
  #           smooth.sd - measure of the sd of noise used in smooth bootstraps
  #           dist.func - character name of the function producing our deviates
  #
  # output  : NULL if no checks fail, otherwise, the error message relevant to 
  #           the first check which failed
  
  if (class(data)!='numeric' & class(data)!='integer'){
    stop('input data must be numeric')}
  
  if ((method=='parametric'| method=='par.fit') &   # if parametric: need a dist
      !(is.character(dist.func) &                   # func name which is a func
        is.function(try(match.fun(dist.func),silent=T))) ){
    stop('when method is parametric or par.fit, dist.func must be provided')
  }
  
  if (method=='smooth' & (class(smooth.sd)!='numeric' | smooth.sd<0 |
                          length(smooth.sd)>1 )){
    stop('When method = \'smooth\', smooth.sd must be a positive scalar')
  }
  
  if (n%%1!=0 | n<2) stop('n must be a positive integer greater than 1')
  
  if (alpha<0 | alpha>1) stop('alpha must be between 0 and 1')
  
  if (!is.function(func)) stop('invalid function supplied as func argument')
  
  if (!method %in% c('percentile','BCa','parametric','smooth','par.fit')){
    stop('invalid method')}
  
}

bootstrap <- function(data, n=999, alpha = 0.05, func = mean,
                      method = 'percentile', smooth.sd = 0.2, 
                      dist.func = NULL, check.inputs=T, ...){
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
  #           method     - 'percentile', 'BCa', 'parametric', 'smooth','par.fit'
  #                        Specifies the bootstrap method to be used. When
  #                        'parametric' is chosen, the percentile method is used
  #                        to calculate the confidence interval using the
  #                        bootstrap samples, and a function from which to
  #                        sample the data must be specified. All remaining 
  #                        options produce non-parametric bootstraps. Option
  #                        'smooth' adds a normal noise centred at 0. The
  #                        chosen standard deviation is a fraction of the sample
  #                        sd. This is set using the parameter 'smooth.sd'.
  #                        'par.fit' is a parametric percentile bootstrap, but
  #                        it estimates the true parameters of the distribution
  #                        from the data, rather than using the known values.
  #                        Only works with dist.func = rnorm, rpois or rgamma.
  #           smooth.sd  - Multiplier for the sample standard deviation. When 
  #                        method = 'smooth', a normal noise id added to 
  #                        each bootstrap resample. It has mean 0 and standard
  #                        deviation smooth.sd * sd(data).
  #           dist.func  - function to sample the data from when parametric is 
  #                        set to TRUE. It is assumed that the first argument
  #                        in any call to dist.func is the number of random
  #                        deviates to be produced, as is the convention with
  #                        rnorm, runif, rpois, rgamma, etc. must be the
  #                        character name of the function.
  #           check.inputs - Logical, if TRUE, all inputs are type checked
  #           ...        - extra optional parameters to be passed to dist.func
  #          
  # output  : named vector containing the lower and upper bounds of the interval
  
  # to allow the user to pass in the func or func name:
  func <- try(match.fun(func), silent=T)
  
  #### Input checks:
    
  if (check.inputs){
    # if dist.func has not been set, call the checking function with NULL
    if (is.null(dist.func)) bootstrap.type.checks(data, n, alpha, func, method,
                                                  smooth.sd, NULL)
    
    # otherwise, pass in dist.func itself. If we tried to do this when
    # dist.func is.null(), we would get an error, hence this unsatisfactory
    # way of dealing with the problem:
    else {bootstrap.type.checks(data, n, alpha, func, method, smooth.sd,
                                dist.func)}
  }
  
  ### End of input-checks
  
  if (is.null(dist.func)==FALSE & (method=='parametric'| method=='par.fit')){
    dist.func.name <- dist.func       # for method=par.fit, we need to know
    dist.func <- match.fun(dist.func) # the name of dist.func
  }
  
  # If the user has set smooth.sd to 0 and wants a smooth bootstrap, we obtain
  # the same result more efficiently by simply providing them with a percentile
  # bootstrap:
  if (smooth.sd==0 & method=='smooth') method <- 'percentile'
  
  if (method!='parametric' & method!='par.fit'){
    
    # generate the random samples with replacement:
    samples <- replicate(n, sample(x=data, size=length(data),replace=T)) 
    
    if (method=='smooth'){ # add noise to the data for a smooth bootstrap
      noise <- replicate(n, rnorm(length(data), sd = sd(data)*smooth.sd))
      samples <- samples + noise
    }
    
  }
  
  else{ # parametric resamples
    
    if (method=='par.fit'){
      
      # get MLEs for the parameters of the distribution:
      dist.func.args <- switch(dist.func.name,
                               
                               # if the data are normal, this is easy using
                               # sample statistics
                               rnorm=list(n=length(data),mean=mean(data),
                                          sd=sd(data)),
                               
                               # same goes for poisson data
                               rpois=list(n=length(data),lambda=mean(data)),
                               
                               # if the data are gamma, we need to call an MLE
                               # function to get the estimates, and extract
                               # initial guesses from the ones passed by the
                               # user using ... arguments. 
                               rgamma=as.list(c(n=length(data),gammaMLE(
                                 log(c(list(...)$rate,list(...)$shape)),data))))
      
      # Note, we add in n=length(data) to the list of estimated parameters so
      # that the use of do.call() becomes possible
    }
    
    # If method=='parametric' we simply use the true parameter values as the 
    # ones we pass to dist.func:
    else {dist.func.args <- as.list(c(n=length(data),list(...)))}
    
    samples <- matrix(nrow=length(data),ncol=n)
    
    # the replicate function worked badly with functions like rpois and rt, 
    # so a less efficient method has to be used to produce parametric samples:
    for (i in 1:n){samples[,i] <- do.call(dist.func,dist.func.args)}
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
        
        dist.function <- match.fun(dist.func)
        
        for (i in 1:simulations){
          dataset <- dist.function(sample.n.setting, ...) # get the O.G. sample
          
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
  
  output <- array(dim = c(dims[1],dims[2],dims[3],3))
  dimnames(output) <- dimnames(simulation.output.object)
  for (i in 1:dims[1]){     # With nested for loops, go through the 
    for(j in 1:dims[2]){    # simulated bootstrap intervals and calculate
      for(k in 1:dims[3]){  # the summary statistics of interest:
        
        boot.ints <- simulation.output.object[i,j,k,]  # extract intervals
        
        coverage <- get.coverage(boot.ints,true.value) # calculate the
        length <- get.length(boot.ints)                # statistics
        failure.tend <- get.coverage(boot.ints,true.value,failure.t=T)
         
        summaries <- c(coverage,length,failure.tend)   # add them to the output
        names(summaries) <- c('coverage','length','failure tendency') #  object
        output[i,j,k,] <- summaries                    # with appropriate names
      }
    }
  }
  
  class(output) <- 'simulation.summary.object'
  return(output)
}

get.coverage <- function(bootstrap.results, true.value, failure.t=FALSE){
  # purpose : returns the observed coverage, given a vector which contains
  #           a sequence of confidence intervals
  #
  # input   : bootstrap.results - a vector containing bootstrap intervals in the
  #                               format c(lower1, upper1, lower2, upper2, etc.)
  #           true.value        - the true value of the statistic of interest. 
  #                               Allows for the calculation of the coverage
  #           failure.t         - failure tendency. Allows the function to
  #                               return the failure tendency rather than the 
  #                               coverage. Failure tendency is a measure (from
  #                               0 to 1), of the proportion of the time the
  #                               true value of the statistic was to the left of
  #                               the confidence interval
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
  
  if(!failure.t){return(sum(in.interval)/(n/2))} # return the observed coverage
  
  else{
    not.in <- as.logical(1-in.interval)
    failed.lowers <- lowers[not.in]                   # return the observed 
    return(sum(true.value<failed.lowers)/sum(not.in)) # failure tendency
  }
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
                                           statistic='coverage',...){
  # purpose : plots the statistic of interest for a set of simulation
  #           bootstrap confidence intervals, for all levels of 'factor'. Fixes
  #           the other setting values at their highest setting i.e. uses the 
  #           the largest sample size and bootstrap resamples available.
  #
  # inputs  : simulation.summary.object - array of summary statistics for 
  #                                       simulation intervals.
  #           statistic                 - summary statistic of interest, 
  #                                       'coverage', 'length','failure
  #                                       tendency'
  #           ...                       - extra optional parameters to be 
  #                                       passed to matplot
  # output  : None, produces a plot. 
  
  if (class(simulation.summary.object)!='simulation.summary.object'){
    stop('invalid input type')}
  
  if ( !(statistic %in% c('coverage','length','failure tendency')) ){
    stop('invalid choice of statistic')}
  
  # fetch summary statistic index:
  stat.ind <- switch(statistic,'coverage'=1, 'length'=2, 'failure tendency'=3) 
  
  # get dimensions of summary object:
  dims <- dim(simulation.summary.object)
  dims.not.stats <- dim(simulation.summary.object[,,,1])
  
  msg1 <- 'Can only plot summaries when sample.n, bootstrap.n and'
  msg2 <- 'method are all vectors'
  if (any(dims.not.stats<2)) stop(paste(msg1, msg2))
  
  for (plot.num in c(1,2)){
  # generate sample size plot first, then bootstrap plot
      
    # extract x axis values from the simulation.summary.object dimnames:
    x <- as.numeric(gsub('[^0-9]','',
                         dimnames(simulation.summary.object)[[plot.num]]))
      
    # extract statistic values and average over index not being plotted:
    y <- apply(simulation.summary.object[,,,stat.ind], c(plot.num,3), mean)
      
    xlab = c('sample size','bootstrap resamples')[plot.num]
  
    # Draw the plot:
    method.names <- gsub('boot.method: ','',
                         dimnames(simulation.summary.object)[[3]])

    matplot(x,y,ylab=statistic,xlab=xlab,type='l',col=seq(1,dims[3]),...)
  
    legend('topright',method.names,lty=1,col=seq(1,dims[3]),bty='n',cex=.75)
  }
}

sim.plot.3D <- function(simulation.summary.object, statistic, method,hist=F,
                        ...){
  # purpose : produces 3D plots of summary statistics for simulated bootstrap 
  #           intervals. Uses methods available by various packages and produces 
  #           2 different types of 3D plot.
  #
  # inputs  : simulation.summary.object - the multi-dimensional array containing
  #                                       the summary statistics for each level
  #                                       of simulation setting.
  #           statistic                 - the character name of the statistic
  #                                       to be plotted
  #           method                    - the integer index of the method to be
  #                                       plotted. Which values are valid
  #                                       depends on the shape of the object 
  #                                       passed to calculate.summaries
  #           hist                      - logical parameter. If TRUE, produces a 
  #                                       3D histogram instead of a 3D
  #                                       3D perspective plot
  #           ...                       - extra optional parameters to be passed
  #                                       to the scatter 3D function.
  #
  # ouput   :  list containing:
  #            x - the sample size values used for the persp3D plot
  #            y - the bootstrap resample values used for the persp3D plot
  #            z - the matrix of statistic values corresponding to the z
  #                values of x and y for the persp3D plot
  #
  # note    : The output is likely of no use to the user, unless they choose
  #           to obtain a plot using this data and a different 3D plotting 
  #           function. The purpose of the output is primarily for debugging
  #           and ensuring the data look as expected. 
  
  ### input checks:
  
  if (class(simulation.summary.object)!='simulation.summary.object'){
    stop('invalid input type')}
  
  if ( !(statistic %in% c('coverage','length','failure tendency')) ){
    stop('invalid choice of statistic')}
  
  if (method<1 | method%%1!=0 | method>dim(simulation.summary.object)[3]){
    stop('invalid choice of method')}
  
  if (any(dim(simulation.summary.object[,,,1])<2)){
    stop('all simulation settings must be vectors to plot in 3D')
  }
  
  ### end of input checks ###
  
  
  # fetch summary statistic index:
  stat.ind <- switch(statistic, 'coverage'=1, 'length'=2, 'failure tendency'=3)
  
  Dnames <- dimnames(simulation.summary.object)               # extract method
  method.name <- gsub('boot.method: ','',Dnames[[3]][method]) # name
  
  ### Format the data for the call to scatter3D and produce a 3D scatter:
  M <- melt(simulation.summary.object[,,method,stat.ind])
  
  x <- as.numeric(gsub('[^0-9]','',M$Var1)) # extract x, y, and z coordinate
  y <- as.numeric(gsub('[^0-9]','',M$Var2)) # values from the melted object
  z <- M$value                              # and produce our plot
  
  scatter3D(x, y, z, main=method.name,xlab='sample size',
            ylab='bootstrap resamples', zlab=statistic,...)
  
  ### Format the data for the call to persp3D and draw the surface:
  x <- as.numeric(gsub('sample.n: ','',Dnames[[1]])) # extract the numeric 
  y <- as.numeric(gsub('boot.n: ','',Dnames[[2]]))   # values of the xs and ys, 
  z <- simulation.summary.object[,,method,stat.ind]  # and the matching z matrix
  
  ifelse(hist, func3D <- hist3D, func3D <- persp3D)
  
  func3D(x,y,z,xlab='sample size',ylab='bootstrap resamples', zlab=statistic,
          main=method.name)
  
  invisible(list(x,y,z)) # to avoid a potentially large matrix from printing
}

gamma.neg.log.lik <- function(par, x){
  # purpose : evaluates the negative log likelihood of a gamma distribution
  #           given parameter guesses on the real line, and data x
  #
  # inputs  : par - parameter estimates for rate and shape as a vector on the
  #                 real line. A log link is applied to transform these values
  #                 to positive ones. 
  #
  # output  : numeric scalar, the negative log likelihood evaluated at x and 
  #           the transformed par
  
  par <- exp(par) # log links to keep alpha and beta positive
  alpha <- par[1] ; beta <- par[2]
  loglik <- dgamma(x, rate=alpha, shape=beta) %>% log %>% sum
  return(-loglik)
}

gammaMLE <- function(par,x,...){
  # purpose : Maximum likelihood estimation of parameters for a gamma
  #           distribution, given observations and initial guesses on the 
  #           real line
  #
  # inputs  : par - values such that exp(par) gives the initial estimates of 
  #                 the rate and shape parameters of the gamma distribution, 
  #                 respectively
  #           x   - vector of observations from the gamma process in question
  #           ... - extra optional parameters to be passed to optim
  #
  # output  : the estimated parameters as a list
  ests <- exp(optim(par, gamma.neg.log.lik, x=x,...)$par)
  return(list(rate=ests[1],shape=ests[2]))
}