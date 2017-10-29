# I confirm that the attached is my own work, except where clearly indicated
# in the text. 

non.parametric.sample <- function(data, n){
  # purpose : produces n random samples of size 'size' from the supplied data
  #
  # inputs  : data - numeric vector of univariate observations
  #           n    - positive integer number of samples to be drawn
  #
  # output  : size*n dimension matrix. Each column is a generated sample.
  
  if (class(data)!='numeric') stop('input data must be numeric')
  if (n%%1!=0 | n<1) stop('n must be a positive integer')
  
  return(replicate(n, sample(x=data, size=length(data),replace=T)))
}

np.p.bootstrap <- function(data, n=999, alpha = 0.05, func = mean,
                           method = 'percentile'){
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
  #          
  # output  : named vector containing the lower and upper bounds of the interval
  
  # Input checks and setting default values: 
  if (class(data)!='numeric') stop('input data must be numeric')
  size <- length(data) 
  if (n%%1!=0 | n<1) stop('n must be a positive integer')
  if (alpha<0 | alpha>1) stop('alpha must be between 0 and 1')
  func <- match.fun(func) # to allow the user to pass in the func or func name
  if (!is.function(func)) stop('invalid function supplied as func argument')
  if (!method %in% c('percentile','bca')) stop('invalid method')
  
  samples <- non.parametric.sample(data, n)       # generate the random samples
  samples <- cbind(samples,data)                  # add in the observed data
  samples <- apply(samples,2,func)                # calculate statistics
  
  lower <- alpha/2      # percentile method
  upper <- 1 - alpha/2  # intervals
  
  if (method=='bca'){
    # calculate BCA intervals
  }
  
  CI <- quantile(samples, probs=c(lower, upper))
  return(CI)
  
}